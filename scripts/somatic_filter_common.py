"""
Shared pieces of the 3-stage somatic call filter (stdlib + PyYAML only):

    basic_filter.py       01 raw -> 02 preprocess (+ on target) -> 03 MANE -> 04-07 sequencing  (vep.csv -> .basic_filtered.csv)
    recurrent_filter.py   08 recurrent artifacts -> 08b indel alt depth   (.basic_filtered.csv -> .recurrent_filtered.csv)
    annotation_filter.py  09a LoF4 -> 09b gnomAD -> 09c germline tag      (.recurrent_filtered.csv -> .annotation_filtered.csv)

All three read config/somatic_filter.yaml (--config + --caller). Each stage writes its own
PREFIX.<stage>.qc.tsv.gz keyed by `row_id` (input row number, assigned in stage 1 and carried through), w/
last_step = index into STEPS below (global across stages), so the R reports can stitch the stages back
into one row-per-raw-call table (scripts/somatic_filter_common.R).
"""
import os
import re
import csv
import time
import argparse

STEPS=['01 raw','02 preprocessed','03 MANE Select','04 orientation','05 depth','06 alt depth','07 VAF',
       '08 recurrent removed','08b indel alt depth','09a LoF4 removed','09b gnomAD removed','09c germline removed']
# last_step values -- index of the last step a row survived
I_SEQ=STEPS.index('07 VAF')                   # passed basic -> enters recurrent_filter.py
I_REC08=STEPS.index('08 recurrent removed')   # passed 08, failed 08b
I_REC=STEPS.index('08b indel alt depth')      # passed 08 + 08b -> enters annotation_filter.py
I_LOF4=STEPS.index('09a LoF4 removed')
I_GNOMAD=STEPS.index('09b gnomAD removed')
I_KEPT=STEPS.index('09c germline removed')    # kept

NUM_COLS=['Sample.AltFrac','Sample.Depth','Sample.AltDepth']
VARIANT_COLS=['Gene','Chr','Start','REF','ALT','Variant.Consequence','Variant.LoF_level','HGVSp','Existing.variation']
QC_COLS=['Sample.ID','Chr','Start','REF','ALT','Gene','FILTER','Variant.LoF_level','Variant.Consequence','Variant.Class']+NUM_COLS

# output file names, PREFIX = input vep csv minus .csv
STAGE_SUFFIX={'basic':'basic_filtered','recurrent':'recurrent_filtered','annotation':'annotation_filtered'}

CSV_OPTS=dict(delimiter=',',restval='.',extrasaction='ignore',quoting=csv.QUOTE_NONNUMERIC,dialect='excel')

# every key config/somatic_filter.yaml may hold (typo guard) -> the stage that uses it
CONFIG_KEYS={'input':'all','output_prefix':'all','label':'report',
             'genes':'basic','targets_bed':'basic','skip_orientation':'basic','min_ad':'basic',
             'min_dp':'basic','vaf_min':'basic','vaf_max':'basic',
             'remove_pct':'recurrent','flag_pct':'recurrent','min_ad_indel':'recurrent',
             'protected':'recurrent','protected_genes':'recurrent','protected_min_vaf':'recurrent','chip_whitelist':'recurrent',
             'gnomad':'annotation','gnomad_exclude_pops':'annotation','gnomad_max_af':'annotation'}

# `gnomad:` config setting -> AF columns to take the max over (vep_vcf_parser.py, gnomAD v4.1 VEP cache).
# 'AF' = the dataset's overall AF (all populations pooled); the rest are per-population AFs.
_GNOMAD_POPS={'gnomADe':['AF','AFR','AMR','ASJ','EAS','FIN','MID','NFE','REMAINING','SAS'],
              'gnomADg':['AF','AFR','AMI','AMR','ASJ','EAS','FIN','MID','NFE','REMAINING','SAS']}
GNOMAD_MODES=('exome','genome','both','max_af')


def gnomad_cols_for(mode,exclude_pops=()):
    """exome / genome / both -> overall AF + per-population AF cols, minus exclude_pops (e.g. bottlenecked
    FIN / ASJ / AMI / REMAINING); max_af -> VEP's MAX_AF col (gnomAD + 1000 Genomes, can't exclude pops)."""
    if mode not in GNOMAD_MODES:
        raise SystemExit(f'gnomad: {mode!r} -- must be one of {list(GNOMAD_MODES)}')
    if mode=='max_af':
        return ['MAX_AF']
    exclude={p.strip().upper() for p in exclude_pops if p.strip()}
    unknown=exclude-{p for pops in _GNOMAD_POPS.values() for p in pops}-{'AF'}
    if unknown:
        raise SystemExit(f'gnomad_exclude_pops: unknown population(s) {sorted(unknown)}')
    if 'AF' in exclude:
        raise SystemExit('gnomad_exclude_pops: AF is the overall AF, not a population')
    srcs={'exome':['gnomADe'],'genome':['gnomADg'],'both':['gnomADe','gnomADg']}[mode]
    return [f'{s}.{p}' for s in srcs for p in _GNOMAD_POPS[s] if p not in exclude]


def to_num(v):
    # '.', NA and -1 all mean missing
    try:
        x=float(v)
    except (TypeError,ValueError):
        return None
    return None if x!=x or x==-1 else x


def locus_id(row):
    return f"{row['Chr']}_{row['Start']}_{row['REF']}_{row['ALT']}"


def gnomad_max(row,cols):
    """Max AF over cols ('&'-joined values -> their max); None = absent from gnomAD."""
    vals=[to_num(x) for c in cols for x in row[c].split('&')]
    vals=[x for x in vals if x is not None]
    return max(vals) if vals else None


def split_cols(v):
    return [c.strip() for c in (v if isinstance(v,list) else str(v).split(',')) if c.strip()]


def median(xs):
    xs=sorted(x for x in xs if x is not None)
    if not xs:
        return '.'
    m=len(xs)//2
    return xs[m] if len(xs)%2 else (xs[m-1]+xs[m])/2


def write_tsv(path,header,rows):
    with open(path,'w',newline='') as fh:
        w=csv.writer(fh,delimiter='\t',lineterminator='\n')
        w.writerow(header)
        w.writerows(rows)


def read_params(report_dir):
    """A stage's params.tsv -> {param: value} ({} if missing)."""
    path=f'{report_dir}/params.tsv'
    if not os.path.exists(path):
        return {}
    with open(path) as fh:
        return {r['param']:r['value'] for r in csv.DictReader(fh,delimiter='\t')}


def new_run_id():
    return time.strftime('%Y%m%d-%H%M%S')+f'-{os.getpid()}'


def open_rows(path,required):
    """csv.DictReader over path after checking required columns; returns (fieldnames, reader, fh)."""
    fh=open(path,newline='')
    reader=csv.DictReader(fh)
    missing=[c for c in required if c not in (reader.fieldnames or [])]
    if missing:
        fh.close()
        raise SystemExit(f'{path} is missing columns: {missing}')
    return reader.fieldnames,reader,fh


# ---- config -----------------------------------------------------------------------------------------
def load_config(path,caller):
    """config/somatic_filter.yaml -> {key: value}: `defaults:` overlaid by `callers.<caller>:`."""
    try:
        import yaml
    except ImportError:
        raise SystemExit('--config needs PyYAML in this env (conda install -n vep_parser pyyaml)')
    with open(path) as fh:
        cfg=yaml.safe_load(fh) or {}
    merged=dict(cfg.get('defaults') or {})
    callers=cfg.get('callers') or {}
    if caller not in callers:
        raise SystemExit(f'--caller {caller} not in {path} callers: {sorted(callers)}')
    merged.update(callers[caller] or {})
    unknown=sorted(set(merged)-set(CONFIG_KEYS))
    if unknown:
        raise SystemExit(f'{path}: unknown key(s) {unknown}; allowed: {sorted(CONFIG_KEYS)}')
    if merged.get('protected')=='':
        merged['protected']=None
    if 'gnomad_exclude_pops' in merged:
        merged['gnomad_exclude_pops']=','.join(split_cols(merged['gnomad_exclude_pops'] or []))
    if 'protected_genes' in merged:
        merged['protected_genes']=','.join(split_cols(merged['protected_genes'] or []))
    return merged


def stage_parser(description):
    """argparse w/ the args every stage shares; returns (parser, pre-parsed --config/--caller)."""
    pre=argparse.ArgumentParser(add_help=False)
    pre.add_argument('--config',help='config/somatic_filter.yaml')
    pre.add_argument('--caller',help='key under `callers:` in the config')
    p=argparse.ArgumentParser(description=description,parents=[pre])
    p.add_argument('-i','--input',help='Merged vep csv (annotate.smk merge_samples) -- or config `input`. '
                                       'Every stage takes this same file; it only sets PREFIX.')
    p.add_argument('-o','--output_prefix',help='PREFIX for every output (default: input minus .csv).')
    return p,pre


def finish_args(p,pre,argv):
    """Config -> argparse defaults (explicit CLI flags still win), then fill PREFIX."""
    known,_=pre.parse_known_args(argv)
    if known.caller and not known.config:
        p.error('--caller needs --config')
    if known.config:
        if not known.caller:
            p.error('--config needs --caller')
        cfg=load_config(known.config,known.caller)
        # only keys this stage has an option for (the rest belong to other stages)
        dests={a.dest for a in p._actions}
        p.set_defaults(**{k:v for k,v in cfg.items() if k in dests})
    args=p.parse_args(argv)
    if not args.input:
        p.error('-i / config `input` is required')
    if not args.output_prefix:
        args.output_prefix=args.input[:-4] if args.input.endswith('.csv') else args.input
    return args


def stage_paths(prefix,stage):
    s=STAGE_SUFFIX[stage]
    return {'csv':f'{prefix}.{s}.csv','qc':f'{prefix}.{s}.qc.tsv.gz','report':f'{prefix}.{s}_report'}


def check_upstream(prefix,stage,upstream):
    """Upstream stage's csv + params must exist; returns its run_id (recorded in this stage's params.tsv
    so the reports can tell when a stage was rerun without rerunning the ones after it)."""
    up=stage_paths(prefix,upstream)
    if not os.path.exists(up['csv']):
        raise SystemExit(f'{up["csv"]} not found -- run the {upstream} stage first')
    run_id=read_params(up['report']).get('run_id')
    if not run_id:
        raise SystemExit(f'{up["report"]}/params.tsv has no run_id -- rerun the {upstream} stage')
    return up,run_id


# ---- protected variants (recurrent_filter.py) --------------------------------------------------------
AA3={'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C','Gln':'Q','Glu':'E','Gly':'G','His':'H','Ile':'I',
     'Leu':'L','Lys':'K','Met':'M','Phe':'F','Pro':'P','Ser':'S','Thr':'T','Trp':'W','Tyr':'Y','Val':'V',
     'Sec':'U','Ter':'*'}
MISSENSE_RE=re.compile(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})')        # p.His1047Arg, p.Arg123Ter
FS_RE=re.compile(r'p\.([A-Z][a-z]{2})(\d+)(?:[A-Z][a-z]{2})?fs')          # p.Arg806ThrfsTer14


def protein_keys(hgvsp):
    """MANE HGVSp -> keys to look up in the protected set: exact ('H1047R' / 'R603*' / 'N1068fs')
    and, for missense only, the residue ('H1047', matched by wildcard rows like KRAS G12)."""
    m=FS_RE.match(hgvsp)
    if m and m.group(1) in AA3:
        return [f'{AA3[m.group(1)]}{m.group(2)}fs']
    m=MISSENSE_RE.fullmatch(hgvsp)
    if m and m.group(1) in AA3 and m.group(3) in AA3:
        ref,pos,alt=AA3[m.group(1)],m.group(2),AA3[m.group(3)]
        return [f'{ref}{pos}{alt}'] if alt=='*' else [f'{ref}{pos}{alt}',f'{ref}{pos}']
    return []


def load_protected(path):
    """gene<TAB>protein_change<TAB>source, '#' lines skipped -> {(gene, protein_change)}."""
    if not path:
        return set()
    with open(path) as fh:
        lines=[l for l in fh if l.strip() and not l.startswith('#')]
    rows=list(csv.DictReader(lines,delimiter='\t'))
    if rows and not {'gene','protein_change'}<=set(rows[0]):
        raise SystemExit(f'{path}: needs gene + protein_change columns, got {list(rows[0])}')
    return {(r['gene'].strip(),r['protein_change'].strip()) for r in rows}


def is_protected(row,protected,protected_genes=frozenset()):
    """Whole gene in protected_genes, or Gene + protein change in the protected variant list."""
    if row['Gene'] in protected_genes:
        return True
    return bool(protected) and any((row['Gene'],k) in protected for k in protein_keys(row['HGVSp']))


def protect_state(listed,ad,dp,floor):
    """'TRUE' = protection applies, 'BELOW_FLOOR' = on the list but AltDepth/Depth < floor, 'FALSE' = not listed."""
    if not listed:
        return 'FALSE'
    return 'TRUE' if ad is not None and dp and ad/dp>=floor else 'BELOW_FLOOR'
