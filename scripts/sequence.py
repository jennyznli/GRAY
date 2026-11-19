"""
Sequence-level filters + flags on one caller's cohort VEP csv. 

    python scripts/sequence.py 
        -i mutect2/all_samples.mutect2.no_pon.vep.csv \
        -o filter/mutect2_no_pon/sequence.csv 
        --removed filter/mutect2_no_pon/sequence_removed.csv \
        --stats reports/filter/sequence_stats.mutect2_no_pon.tsv \
        --genes input/IDT-Target-HYB-PreCancer.genes.txt 
        --targets-bed <slop50 merged bed>

Removed (-> --removed, + removed_at = the step it failed):
  01  Gene in --genes, Bio.type == protein_coding, AltFrac / Depth / AltDepth present, on --targets-bed
  02  Sample.Depth >= --min-dp
  03  REF and ALT <= --max-len bp
  04  Sample.F1R2 alt >= 1 and Sample.F2R1 alt >= 1 (Mutect2 only)

Flags 
  pass_minad   Sample.AltDepth >= --min-ad-snv, or >= --min-ad-indel when len(REF) != len(ALT)
  pass_strand  Sample.AltFwd / (AltFwd + AltRev) within [--strand-min, 1 - --strand-min]; NA w/o them
"""
import csv
import bisect
import argparse
from collections import Counter,defaultdict

STEPS=['s01_target_coding','s02_dp','s03_indel_len','s04_alt_orientation']
VALUE_COLS=('Sample.AltFrac','Sample.Depth','Sample.AltDepth')


def num(v):
    """float; None for '', '.', NA, -1."""
    try:
        x=float(v)
    except (TypeError,ValueError):
        return None
    return None if x!=x or x==-1 else x


def pair(v):
    """'a|b' -> (a, b); None if absent."""
    try:
        a,b=(float(c) for c in v.split('|'))
    except (AttributeError,ValueError):
        return None
    return a,b


def load_bed(path):
    """BED -> {chrom: (starts, ends)}, merged, 0-based half-open."""
    iv=defaultdict(list)
    with open(path) as fh:
        for l in fh:
            f=l.split('\t')
            if len(f)>=3 and f[1].isdigit():
                iv[f[0]].append((int(f[1]),int(f[2])))
    if not iv:
        raise SystemExit(f'{path}: no tab-separated BED intervals')
    bed={}
    for chrom,ivs in iv.items():
        merged=[]
        for s,e in sorted(ivs):
            if merged and s<=merged[-1][1]:
                merged[-1][1]=max(merged[-1][1],e)
            else:
                merged.append([s,e])
        bed[chrom]=([s for s,_ in merged],[e for _,e in merged])
    return bed


def on_target(x,bed):
    """Start .. Start + len(REF) - 1 overlaps a BED interval."""
    if bed is None:
        return True
    if x['Chr'] not in bed:
        return False
    starts,ends=bed[x['Chr']]
    start=int(x['Start'])
    i=bisect.bisect_left(starts,start+max(len(x['REF']),1)-1)-1
    return i>=0 and ends[i]>start-1


def main(argv=None):
    p=argparse.ArgumentParser()
    p.add_argument('-i','--input',required=True,help='cohort vep csv')
    p.add_argument('-o','--output',required=True)
    p.add_argument('--removed',required=True,help='rows that fail a step')
    p.add_argument('--stats',required=True)
    p.add_argument('--genes',default='',help='one gene per line; "" = no restriction')
    p.add_argument('--targets-bed',default='',help='"" = no restriction')
    p.add_argument('--min-dp',type=float,default=100)
    p.add_argument('--max-len',type=int,default=100)
    p.add_argument('--min-ad-snv',type=float,default=5)
    p.add_argument('--min-ad-indel',type=float,default=10)
    p.add_argument('--strand-min',type=float,default=0.1)
    a=p.parse_args(argv)

    genes=None
    if a.genes:
        with open(a.genes) as fh:
            genes={l.strip() for l in fh if l.strip() and not l.startswith('#')}
    bed=load_bed(a.targets_bed) if a.targets_bed else None

    n=Counter()
    with open(a.input,newline='') as fi, open(a.output,'w',newline='') as fo, open(a.removed,'w',newline='') as fr:
        r=csv.DictReader(fi)
        w=csv.DictWriter(fo,fieldnames=r.fieldnames+['pass_minad','pass_strand'])
        w_rm=csv.DictWriter(fr,fieldnames=r.fieldnames+['removed_at'])
        w.writeheader()
        w_rm.writeheader()
        for x in r:
            n['input']+=1
            dp,ad=num(x['Sample.Depth']),num(x['Sample.AltDepth'])
            f1,f2=pair(x.get('Sample.F1R2')),pair(x.get('Sample.F2R1'))
            # genes is in list (if given) and protein coding 
            ok=[(genes is None or x['Gene'] in genes) and x['Bio.type']=='protein_coding'
                and all(num(x[c]) is not None for c in VALUE_COLS) and on_target(x,bed),
                dp is not None and dp>=a.min_dp,
                # max length of REF or ALT is less than 100  
                max(len(x['REF']),len(x['ALT']))<=a.max_len,
                # Mutect2 - both alt reads are non zero 
                (f1 is None and f2 is None) or (f1 is not None and f2 is not None and f1[1]>=1 and f2[1]>=1)]
            cum=True
            for step,s_ok in zip(STEPS,ok):
                cum=cum and s_ok
                n[step]+=cum
            if not cum:
                x['removed_at']=STEPS[ok.index(False)]
                w_rm.writerow(x)
                continue
            # count rows with no orientation info
            n['no_orientation']+=f1 is None and f2 is None

            # flag indels/snv that don't pass min ad 
            indel=len(x['REF'])!=len(x['ALT'])
            x['pass_minad']=str(ad>=(a.min_ad_indel if indel else a.min_ad_snv)).upper()
            
            # flag rows that don't pass strand filtering
            s=[num(x.get('Sample.AltFwd')),num(x.get('Sample.AltRev'))]
            x['pass_strand']='NA' if None in s or sum(s)==0 else \
                str(a.strand_min<=s[0]/sum(s)<=1-a.strand_min).upper()
                
            n['pass_minad']+=x['pass_minad']=='TRUE'
            n['pass_strand']+=x['pass_strand']=='TRUE'
            n['strand_NA']+=x['pass_strand']=='NA'
            w.writerow(x)

    order=['input']+STEPS+['no_orientation','pass_minad','pass_strand','strand_NA']
    with open(a.stats,'w') as fs:
        fs.write('step\tn\n'+''.join(f'{k}\t{n[k]}\n' for k in order))
    print(f'{a.input}: '+', '.join(f'{k}={n[k]}' for k in order)+f' -> {a.output}')


if __name__=='__main__':
    main()
