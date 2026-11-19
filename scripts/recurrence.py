"""
Stage 2 of 3 (see somatic_filter_common.py): recurrent-artifact filter on basic_filter.py's output.

    python scripts/recurrent_filter.py --config config/somatic_filter.yaml --caller mutect2_no_pon
    # reads PREFIX.basic_filtered.csv (PREFIX from the config `input`, same as every stage)

08 recurrent: per exact variant (Chr:Start:REF:ALT), distinct samples carrying it among the basic-passed rows
  (= every row of the input here). Cutoffs are a fraction of N = samples in this caller's raw input (basic stage
  01 raw n_samples), so callers w/ different sample counts get the same relative cutoff:
  carriers > --remove-pct x N removed; > --flag-pct x N (and not removed) kept w/ flag_recurrent=TRUE.
  e.g. remove_pct 0.5, N = 119 -> > 59.5 -> removed at >= 60 carriers.
Exceptions skip 08: whole genes (--protected-genes, e.g. PIK3CA,TP53) and listed variants (--protected,
Gene + MANE HGVSp) -- only in rows w/ AltDepth/Depth >= --protected-min-vaf (0 = every row; protected = TRUE,
listed but below = BELOW_FLOOR). At a protected variant over the remove cutoff, the >= floor rows are kept w/
flag_recurrent=TRUE and the below-floor rows are removed (carriers still count every sample).
CHIP whitelist (--chip-whitelist DIR, scripts/chip_whitelist.py): rows w/ wl_whitelist or wl_manualreview =
TRUE count as listed too (same VAF floor); the wl_* columns are added to every output row.
The `protected` column written here is what annotation_filter.py uses for the 09b gnomAD exemption.

08b indel alt depth: rows that survive 08 w/ Variant.Class insertion / deletion / indel need AltDepth >=
  --min-ad-indel (SNV / substitution already met basic min_ad). Applied AFTER carriers are counted, so low-AD
  indels still count toward a variant's recurrence. Protection does not exempt rows from 08b.

Two streaming passes (pass 1 counts carriers, pass 2 classifies + writes), so memory stays flat.

Outputs:
  PREFIX.recurrent_filtered.csv          rows passing 08 + n_carriers / flag_recurrent / protected -> annotation_filter.py
  PREFIX.recurrent_removed.csv           rows dropped by 08 (same cols)
  PREFIX.indel_ad_removed.csv            rows dropped by 08b (same cols)
  PREFIX.recurrent_filtered.qc.tsv.gz    every input row: row_id + protected / n_carriers / flag_recurrent / last_step (6-8)
  PREFIX.recurrent_filtered_report/
      params.tsv               thresholds (pct + the absolute carrier cutoffs they give) + run_id + upstream run_id
      funnel.tsv               the 08 + 08b rows of the funnel
      recurrent_variants.tsv   one row per variant, most carriers first
      protected_rescued.tsv    protected variants 08 would otherwise have removed
"""
import os
import csv
import gzip
from collections import Counter,defaultdict

from somatic_filter_common import (STEPS,I_SEQ,I_REC08,I_REC,VARIANT_COLS,CSV_OPTS,to_num,locus_id,median,write_tsv,
                                   read_params,new_run_id,open_rows,stage_parser,finish_args,stage_paths,
                                   check_upstream,load_protected,is_protected,protect_state)
import chip_whitelist

REQUIRED=['row_id','Sample.ID','FILTER','Sample.AltFrac','Sample.Depth','Sample.AltDepth','Variant.Class']+VARIANT_COLS
INDEL_CLASSES={'insertion','deletion','indel'}     # VEP Variant.Class; 'indel' = complex (VarDict)


def get_args(argv):
    p,pre=stage_parser('Stage 2: recurrent-artifact filter (08).')
    p.add_argument('--remove-pct',type=float,default=0.5,help='carriers > this fraction of N samples -> removed.')
    p.add_argument('--flag-pct',type=float,default=None,
                   help='carriers > this fraction of N (and not removed) -> flag_recurrent=TRUE (unset = no flagging).')
    p.add_argument('--min-ad-indel',type=float,default=10,
                   help='08b: indels (Variant.Class insertion / deletion / indel) need AltDepth >= this.')
    p.add_argument('--protected-genes',default='',help='Comma-separated genes exempt from 08 (and 09b) as a whole.')
    p.add_argument('--protected',default=None,help='Variant list (gene / protein_change) exempt from 08 (and 09b).')
    p.add_argument('--protected-min-vaf',type=float,default=0,
                   help='Protection only applies in rows w/ AltDepth/Depth >= this (0 = always).')
    p.add_argument('--chip-whitelist',default=None,
                   help='CHIP whitelist dir (chip_whitelist.py); whitelist + manual-review hits are protected too.')
    args=finish_args(p,pre,argv)
    for name in ('remove_pct','flag_pct'):
        v=getattr(args,name)
        if v is not None and not 0<=v<=1:
            p.error(f'{name} must be a fraction between 0 and 1, got {v}')
    if args.flag_pct is not None and args.flag_pct>args.remove_pct:
        p.error(f'flag_pct ({args.flag_pct}) > remove_pct ({args.remove_pct})')
    return args


def caller_n_samples(report_dir):
    """N = samples in the caller's raw input (basic stage funnel, 01 raw)."""
    with open(f'{report_dir}/funnel.tsv') as fh:
        for r in csv.DictReader(fh,delimiter='\t'):
            if r['step']==STEPS[0]:
                return int(r['n_samples'])
    raise SystemExit(f'{report_dir}/funnel.tsv has no "{STEPS[0]}" row')


def main(argv=None):
    args=get_args(argv)
    up,up_run_id=check_upstream(args.output_prefix,'recurrent','basic')
    out=stage_paths(args.output_prefix,'recurrent')
    os.makedirs(out['report'],exist_ok=True)
    protected=load_protected(args.protected)
    protected_genes=frozenset(c.strip() for c in (args.protected_genes or '').split(',') if c.strip())
    wl=chip_whitelist.Whitelist(args.chip_whitelist) if args.chip_whitelist else None
    required=REQUIRED+([c for c in chip_whitelist.REQUIRED if c not in REQUIRED] if wl else [])
    wl_cols=chip_whitelist.WL_COLS if wl else []
    run_id=new_run_id()
    # fractions -> absolute carrier cutoffs: carriers > pct x N  <=>  carriers > floor(pct x N) (carriers are integers)
    n_samples=caller_n_samples(up['report'])
    remove_threshold=int(args.remove_pct*n_samples+1e-9)      # +1e-9: 0.29*100 = 28.999...
    flag_threshold=remove_threshold if args.flag_pct is None else int(args.flag_pct*n_samples+1e-9)

    def prot_of(row):
        listed=is_protected(row,protected,protected_genes)
        if wl:
            ann=wl.annotate(row)
            row.update(ann)
            listed=listed or ann['wl_whitelist']=='TRUE' or ann['wl_manualreview']=='TRUE'
        return protect_state(listed,to_num(row['Sample.AltDepth']),to_num(row['Sample.Depth']),args.protected_min_vaf)

    # ---- pass 1: carriers per variant
    carriers=defaultdict(set)
    var_info={}
    var_pass=Counter()
    var_vafs=defaultdict(list)
    protected_loci=set()               # variants w/ >= 1 row where protection applies (TRUE)
    prot_rows=Counter()                # variant -> TRUE rows
    fields,reader,fh=open_rows(up['csv'],required)
    with fh:
        for row in reader:
            k=locus_id(row)
            carriers[k].add(row['Sample.ID'])
            if k not in var_info:
                var_info[k]=[row[c] for c in VARIANT_COLS]
            var_pass[k]+=row['FILTER']=='PASS'
            var_vafs[k].append(to_num(row['Sample.AltFrac']))
            if prot_of(row)=='TRUE':
                protected_loci.add(k)
                prot_rows[k]+=1

    n_carriers={k:len(v) for k,v in carriers.items()}
    over_remove={k for k,n in n_carriers.items() if n>remove_threshold}
    removed=over_remove-protected_loci        # every row removed
    rec_rescued=over_remove&protected_loci    # only the protected (>= floor) rows kept
    flagged={k for k,n in n_carriers.items() if flag_threshold<n<=remove_threshold}|rec_rescued

    # ---- pass 2: classify + write
    n_in=n_kept=n_kept_pass=n_removed=n_removed_pass=n_below_floor=n_rec_below_floor=0
    n_indel_removed=n_indel_removed_pass=0
    kept_variants=set()
    kept_samples=set()
    rec_n=rec_pass=0                   # funnel row 08 (passed 08, before 08b)
    rec_variants=set()
    rec_samples=set()
    # n_carriers / sample_freq already come from basic_filter.py (same counts) -> only add what's missing
    out_fields=fields+[c for c in ['n_carriers','flag_recurrent','protected']+wl_cols if c not in fields]
    fields,reader,fh=open_rows(up['csv'],required)
    with fh, open(out['csv'],'w',newline='') as f_kept, \
            open(f'{args.output_prefix}.recurrent_removed.csv','w',newline='') as f_rem, \
            open(f'{args.output_prefix}.indel_ad_removed.csv','w',newline='') as f_ind, \
            gzip.open(out['qc'],'wt',newline='',compresslevel=3) as f_qc:
        w_kept=csv.DictWriter(f_kept,fieldnames=out_fields,**CSV_OPTS)
        w_rem=csv.DictWriter(f_rem,fieldnames=out_fields,**CSV_OPTS)
        w_ind=csv.DictWriter(f_ind,fieldnames=out_fields,**CSV_OPTS)
        w_kept.writeheader()
        w_rem.writeheader()
        w_ind.writeheader()
        w_qc=csv.writer(f_qc,delimiter='\t',lineterminator='\n')
        w_qc.writerow(['row_id','protected','n_carriers','flag_recurrent','last_step'])
        for row in reader:
            n_in+=1
            k=locus_id(row)
            prot=prot_of(row)
            n_below_floor+=prot=='BELOW_FLOOR'
            # removal is per row: at a rescued variant only the protected (>= floor) rows stay
            drop=k in over_remove and prot!='TRUE'
            n_rec_below_floor+=drop and k in rec_rescued
            row['n_carriers']=n_carriers[k]
            row['flag_recurrent']='TRUE' if k in flagged else 'FALSE'
            row['protected']=prot
            is_pass=row['FILTER']=='PASS'
            ad=to_num(row['Sample.AltDepth'])
            drop_indel=not drop and row['Variant.Class'] in INDEL_CLASSES and (ad is None or ad<args.min_ad_indel)
            w_qc.writerow([row['row_id'],prot,n_carriers[k],row['flag_recurrent'],
                           I_SEQ if drop else I_REC08 if drop_indel else I_REC])
            if not drop:
                rec_n+=1
                rec_pass+=is_pass
                rec_variants.add(k)
                rec_samples.add(row['Sample.ID'])
            if drop:
                n_removed+=1
                n_removed_pass+=is_pass
                w_rem.writerow(row)
            elif drop_indel:
                n_indel_removed+=1
                n_indel_removed_pass+=is_pass
                w_ind.writerow(row)
            else:
                n_kept+=1
                n_kept_pass+=is_pass
                kept_variants.add(k)
                kept_samples.add(row['Sample.ID'])
                w_kept.writerow(row)

    # ---- reports
    up_params=read_params(up['report'])
    n_raw=None
    with open(f'{up["report"]}/funnel.tsv') as fh:
        for r in csv.DictReader(fh,delimiter='\t'):
            if r['step']==STEPS[0]:
                n_raw=int(r['n_rows'])
    write_tsv(f'{out["report"]}/params.tsv',['param','value'],
              [['stage','recurrent'],['run_id',run_id],['basic_run_id',up_run_id],
               ['config',args.config or 'none'],['caller',args.caller or 'none'],['input_csv',up['csv']],
               ['n_samples',n_samples],['remove_pct',f'{args.remove_pct:g}'],
               ['flag_pct','none' if args.flag_pct is None else f'{args.flag_pct:g}'],
               ['remove_threshold',remove_threshold],['flag_threshold',flag_threshold],
               ['min_ad_indel',f'{args.min_ad_indel:g}'],
               ['protected_genes',','.join(sorted(protected_genes)) or 'none'],
               ['protected_file',args.protected or 'none'],['n_protected_entries',len(protected)],
               ['protected_min_vaf',f'{args.protected_min_vaf:g}'],
               ['chip_whitelist',', '.join(wl.paths.values()) if wl else 'none']])
    write_tsv(f'{out["report"]}/funnel.tsv',['step','n_rows','pct_of_raw','n_unique_variants','n_PASS_rows','n_samples'],
              [[STEPS[I_REC08],rec_n,round(100*rec_n/n_raw,1) if n_raw else '.',len(rec_variants),rec_pass,
                len(rec_samples)],
               [STEPS[I_REC],n_kept,round(100*n_kept/n_raw,1) if n_raw else '.',len(kept_variants),n_kept_pass,
                len(kept_samples)]])

    def status(k):
        if k in removed:
            return f'removed (> {remove_threshold})'
        if k in rec_rescued:
            return f'protected (> {remove_threshold}, rows >= {args.protected_min_vaf:g} VAF kept)'
        if k in flagged:
            return f'flagged ({flag_threshold+1}-{remove_threshold})'
        return f'kept (<= {flag_threshold})'

    write_tsv(f'{out["report"]}/recurrent_variants.tsv',
              ['locus_id','n_carriers']+VARIANT_COLS+['n_PASS','median_VAF','recurrent_status','protected'],
              sorted(([k,n_carriers[k]]+var_info[k]+[var_pass[k],median(var_vafs[k]),status(k),
                       'TRUE' if k in protected_loci else 'FALSE'] for k in carriers),key=lambda r:-r[1]))
    write_tsv(f'{out["report"]}/protected_rescued.tsv',['locus_id','rescued_from','why','n_rows']+VARIANT_COLS,
              sorted([k,'08 recurrent',f'n_carriers {n_carriers[k]}, kept {prot_rows[k]} rows >= '
                      f'{args.protected_min_vaf:g} VAF',prot_rows[k]]+var_info[k] for k in rec_rescued))

    # ---- stdout summary
    print(f'input: {up["csv"]} (basic run {up_run_id}, {up_params.get("caller","?")})')
    print(f'recurrent: N = {n_samples} samples; remove > {args.remove_pct:g} x N (> {remove_threshold} carriers); '
          f'flag {"off" if args.flag_pct is None else f"> {args.flag_pct:g} x N (> {flag_threshold} carriers)"}')
    print(f'protected: genes {",".join(sorted(protected_genes)) or "none"} + {len(protected)} variant entries from '
          f'{args.protected or "none"}, applied in rows w/ VAF >= {args.protected_min_vaf:g}')
    if wl:
        print(f'CHIP whitelist: {", ".join(wl.paths.values())}')
    print(f'rows in {n_in}  kept {n_kept}  removed {n_removed} ({n_removed_pass} PASS) by 08, '
          f'{n_indel_removed} ({n_indel_removed_pass} PASS) by 08b indel AltDepth < {args.min_ad_indel:g}')
    print(f'variants: {len(removed)} removed, {len(flagged)} flagged ({sum(n_carriers[k] for k in flagged)} rows), '
          f'{len(rec_rescued)} protected kept despite > {remove_threshold} carriers '
          f'({sum(prot_rows[k] for k in rec_rescued)} rows >= floor kept, {n_rec_below_floor} below-floor rows removed)')
    print(f'listed-but-below-floor rows (no protection): {n_below_floor}; protected variants: {len(protected_loci)}')
    print(f'\n{out["csv"]} written ({n_kept} rows); qc + reports next to it')


if __name__=='__main__':
    main()
