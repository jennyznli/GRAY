"""
Merge the callers' sequence.py outputs into one row per call (Sample.ID + variant + Gene).

    python scripts/merge.py
        --mutect2 filter/mutect2_no_pon/sequence.csv \
        --vardict filter/vardict/sequence.csv \
        -o filter/merged.csv
        --removed filter/merged_removed.csv
        --stats reports/filter/merged_stats.tsv

Kept: called by both callers (--both), and at least one caller passes each of (not necessarily the same caller)
  pass_filter  Mutect2: no FILTER tag other than PASS / --mutect2-ignore; VarDict: FILTER == PASS
  pass_minad   sequence.py
  pass_strand  sequence.py
The rest go to --removed (same columns).

Columns
  annotation                                      from Mutect2 when it has the call, else VarDict
  Sample.Depth / AltDepth / AltFrac               mean over the callers that have the call (AltFrac = mean of AltDepth / Depth)
  callers                                         Mutect2;VarDict / Mutect2 / VarDict
  <Caller>.<col>                                  per caller ('.' / NA when the caller lacks the call)
  pass_filter, pass_minad, pass_strand            TRUE in at least one caller
"""
import csv
import argparse
from collections import Counter

# input column -> per-caller name
PER_CALLER={'FILTER':'FILTER','Sample.Depth':'Depth','Sample.AltDepth':'AltDepth','Sample.AltFrac':'AltFrac',
            'Sample.AF':'AF','Sample.F1R2':'F1R2','Sample.F2R1':'F2R1','Sample.AltFwd':'AltFwd',
            'Sample.AltRev':'AltRev','FORMAT':'FORMAT','pass_minad':'pass_minad','pass_strand':'pass_strand'}
SHARED=('Sample.Depth','Sample.AltDepth','Sample.AltFrac')     # also kept under their own name, averaged
FLAGS=('pass_filter','pass_minad','pass_strand')
KEY=('Sample.ID','Chr','Start','REF','ALT','Gene')


def num(v):
    try:
        return float(v)
    except (TypeError,ValueError):
        return None


def mean(vals):
    v=[x for x in vals if x is not None]
    return f'{sum(v)/len(v):.4f}'.rstrip('0').rstrip('.') if v else '.'


def caller_cols(caller):
    names=[n for n in PER_CALLER.values() if not (caller=='VarDict' and n in ('F1R2','F2R1'))]
    return [f'{caller}.{n}' for n in names[:1]+['pass_filter']+names[1:]]


def main(argv=None):
    p=argparse.ArgumentParser()
    p.add_argument('--mutect2',help='sequence.csv')
    p.add_argument('--vardict',help='sequence.csv')
    p.add_argument('-o','--output',required=True)
    p.add_argument('--removed',required=True,help='calls that fail')
    p.add_argument('--stats',required=True)
    p.add_argument('--both',action='store_true',help='keep only calls made by both callers')
    p.add_argument('--mutect2-ignore',default='germline,clustered_events,multiallelic',
                   help='Mutect2 FILTER tags that do not fail a call')
    a=p.parse_args(argv)
    paths={c:f for c,f in (('Mutect2',a.mutect2),('VarDict',a.vardict)) if f}
    if not paths:
        p.error('need --mutect2 and / or --vardict')
    ignore={'PASS'}|set(filter(None,a.mutect2_ignore.split(',')))

    n=Counter()         # (source, column) -> calls
    n_samples={}
    calls={}            # key -> {caller: row}
    base_cols=[]
    for caller,path in paths.items():
        samples=set()
        with open(path,newline='') as fh:
            r=csv.DictReader(fh)
            base_cols+=[c for c in r.fieldnames if c not in base_cols and (c not in PER_CALLER or c in SHARED)]
            for x in r:
                tags=set(x['FILTER'].split(';'))
                x['pass_filter']=str(not tags-ignore if caller=='Mutect2' else tags=={'PASS'}).upper()
                calls.setdefault(tuple(x[k] for k in KEY),{})[caller]=x
                samples.add(x['Sample.ID'])
        n_samples[caller]=len(samples)

    out_cols=base_cols+['callers']+[c for caller in paths for c in caller_cols(caller)]+list(FLAGS)
    with open(a.output,'w',newline='') as fo, open(a.removed,'w',newline='') as fr:
        w=csv.DictWriter(fo,fieldnames=out_cols,extrasaction='ignore')
        w_rm=csv.DictWriter(fr,fieldnames=out_cols,extrasaction='ignore')
        w.writeheader()
        w_rm.writeheader()
        for rows in calls.values():
            src='both' if len(rows)==2 else f'{next(iter(rows))}_only'
            out=dict(next(iter(rows.values())))            # Mutect2's row when present
            out['callers']=';'.join(rows)
            dp=[num(x['Sample.Depth']) for x in rows.values()]
            ad=[num(x['Sample.AltDepth']) for x in rows.values()]
            out['Sample.Depth'],out['Sample.AltDepth']=mean(dp),mean(ad)
            out['Sample.AltFrac']=mean(x/d if x is not None and d else None for x,d in zip(ad,dp))
            for caller in paths:
                x=rows.get(caller,{})
                out[f'{caller}.pass_filter']=x.get('pass_filter','NA')
                for col,name in PER_CALLER.items():
                    out[f'{caller}.{name}']=x.get(col,'NA' if name.startswith('pass_') else '.')
            for f in FLAGS:
                out[f]=str(any(x.get(f)=='TRUE' for x in rows.values())).upper()
            keep=all(out[f]=='TRUE' for f in FLAGS) and (len(rows)==len(paths) or not a.both)
            for s in (src,'all'):
                n[s,'calls']+=1
                n[s,'kept']+=keep
                for f in FLAGS:
                    n[s,f]+=out[f]=='TRUE'
            (w if keep else w_rm).writerow(out)

    # stats: calls per caller source (rows) passing each flag / kept
    cols=['calls',*FLAGS,'kept']
    sources=(['both'] if len(paths)==2 else [])+[f'{c}_only' for c in paths]+['all']
    table='\t'.join(['source']+cols)+'\n'+''.join('\t'.join([s]+[str(n[s,c]) for c in cols])+'\n' for s in sources)
    with open(a.stats,'w') as fs:
        fs.write(table)
    print('samples: '+', '.join(f'{c}={n_samples[c]}' for c in paths))
    print(table+f'-> {a.output}')


if __name__=='__main__':
    main()
