"""
Whitelist flags on the merged calls: CHIP whitelist + tumour hotspots.

    python scripts/whitelist.py
        -i filter/merged.csv \
        -o filter/whitelist.csv
        --protected filter/whitelist_protected.csv
        --stats reports/filter/whitelist_stats.tsv \
        --whitelist-dir whitelist
        --hotspots whitelist/YYYYMMDD_hotspot_rules.tsv

Whitelist dir (latest file per pattern):
  *_chip_gene_list.txt
  *_chip_missense.txt
  *_chip_lof.txt
  *_chip_splice.txt
  *_chip_whitelist_rules.tsv

Added columns (nothing removed):
  wl_whitelist     meets any CHIP whitelist criteria (CHIP genes only)
  wl_mis           missense variant
  wl_lof           loss-of-function variant
  wl_splice        splice variant
  wl_exception     exception
  wl_manualreview  needs manual review
  wl_rule          custom written fields read in from *_chip_whitelist_rules.tsv
  wl_transcript    MANE.Select == the gene list's Accession
  wl_hotspot       matches a --hotspots rule: gene + variant class + protein position (from HGVSp)
  wl_hotspot_rule  the matching rule(s)
  protected        wl_whitelist, wl_manualreview or wl_hotspot
"""
import os
import re
import csv
import glob
import argparse
from collections import Counter

REQUIRED=['Gene','Variant.Consequence','Amino.acids','Protein.position','EXON','INTRON','MANE.Select']

WL_COLS=['wl_whitelist','wl_mis','wl_lof','wl_splice','wl_exception','wl_manualreview','wl_rule','wl_transcript']

# maps VEP consequence term -> whitelist class ('&' consequences are split)
CLASS_OF={'missense_variant':'missense',
        'stop_gained':'nonsense',
        'frameshift_variant':'frameshift',
        'splice_acceptor_variant':'splice',
        'splice_donor_variant':'splice',
        'inframe_deletion':'inframe_deletion',
        'inframe_insertion':'inframe_insertion'}
CLASS_ALIASES={'lof':{'nonsense','frameshift'},'inframe':{'inframe_deletion','inframe_insertion'}}
VALID_CLASSES=set(CLASS_OF.values())|set(CLASS_ALIASES)

# whitelist file reader, stores rows as a list of dicts
def _rows(path):
    with open(path,newline='') as fh:
        # '"#...' = a comment line Excel quoted on save
        lines=[l.replace('\r','') for l in fh if l.strip() and not l.lstrip('"').startswith('#')]
    return [{k.strip():(v or '').strip() for k,v in r.items()} for r in csv.DictReader(lines,delimiter='\t')]

# gene list -> {gene: accession}; a plain one-gene-per-line file (no header) gives accession ''
def _gene_accessions(path):
    with open(path) as fh:
        lines=[l.strip() for l in fh if l.strip() and not l.startswith('#')]
    if lines and lines[0].split('\t')[0]=='Gene':
        return {r['Gene']:r.get('Accession','') for r in _rows(path)}
    return {l.split('\t')[0]:'' for l in lines}

# returns the latest file matching a glob pattern
def _latest(d,pattern):
    hits=sorted(glob.glob(os.path.join(d,pattern)))
    if not hits:
        raise SystemExit(f'{d}: no file matching {pattern}')
    return hits[-1]

# parses protein position ranges as a list of (lo, hi) tuples
def _ranges(s):
    out=[]
    for part in (p.strip() for p in s.split(',') if p.strip()):
        lo,_,hi=part.partition('-')
        out.append((int(lo),int(hi or lo)))
    return out

# returns the first integer in a string, or None if none is found
def _first_int(s):
    m=re.match(r'\s*(\d+)',s or '')
    return int(m.group(1)) if m else None

# returns the first exon or intron number from an EXON / INTRON field
def _region_num(s):
    return _first_int(re.split(r'[|/]',s or '')[0])

# returns the residue span (lo, hi) of a variant: HGVSp (p.Met326IlefsTer19 -> 326, p.Glu746_Ala750del -> 746-750),
# else Protein.position; None if neither has a position
HGVSP_RE=re.compile(r'p\.[A-Z][a-z]{2}(\d+)(?:_[A-Z][a-z]{2}(\d+))?')
def _aa_span(row):
    m=HGVSP_RE.search(row.get('HGVSp') or '')
    pos=[int(x) for x in (m.groups() if m else re.findall(r'\d+',row['Protein.position'].split('/')[0])) if x]
    return (min(pos),max(pos)) if pos else None

# returns the whitelist classes for a VEP consequence string
def var_classes(consequence):
    return {CLASS_OF[t] for t in (consequence or '').split('&') if t in CLASS_OF}

# turns VEP amino acid and protein position fields into a single-residue protein change (compare w whitelist)
def protein_change(amino_acids,protein_position):
    ref,sep,alt=(amino_acids or '').partition('/')
    pos=(protein_position or '').strip()
    if sep and len(ref)==1 and len(alt)==1 and ref!='-' and alt!='-' and pos.isdigit():
        return f'{ref}{pos}{alt}'
    return None

class Whitelist:
    def __init__(self,wl_dir,rules_path=None):
        self.paths={'genes':_latest(wl_dir,'*_chip_gene_list.txt'),
                    'missense':_latest(wl_dir,'*_chip_missense.txt'),
                    'lof':_latest(wl_dir,'*_chip_lof.txt'),
                    'splice':_latest(wl_dir,'*_chip_splice.txt'),
                    'rules':rules_path or _latest(wl_dir,'*_chip_whitelist_rules.tsv')}
        self.accession=_gene_accessions(self.paths['genes'])
        self.missense={(r['Gene'],r['AAChange']) for r in _rows(self.paths['missense'])}
        self.lof_genes={r['Gene'] for r in _rows(self.paths['lof'])}
        self.splice_genes={r['Gene'] for r in _rows(self.paths['splice'])}
        self.rules=[]

        # parse the custom rules file
        for i,r in enumerate(_rows(self.paths['rules']),1):
            classes=set()

            # split the comma-separated list of variant classes
            for c in (x.strip() for x in r['var_class'].split(',') if x.strip()):
                if c not in VALID_CLASSES:
                    raise SystemExit(f'{self.paths["rules"]} row {i}: unknown var_class {c!r}')
                classes|=CLASS_ALIASES.get(c,{c})

            # check the action
            if r['action'] not in ('whitelist','manualreview'):
                raise SystemExit(f'{self.paths["rules"]} row {i}: action must be whitelist / manualreview')
            
            # create a label for the rule and add a restriction if present
            label=f"{r['gene']}_{r['var_class']}"
            restriction=r['exons'] or r['aa_ranges']
            if restriction:
                label+=f"_{restriction}"

            self.rules.append(dict(gene=r['gene'],classes=classes,exons=_ranges(r['exons']),
                                   aa=_ranges(r['aa_ranges']),action=r['action'],
                                   label=label))
                                
        unknown=sorted({r['gene'] for r in self.rules}-set(self.accession))
        if unknown:
            raise SystemExit(f'{self.paths["rules"]}: genes not in the gene list: {unknown}')

        # gene + class pairs a restricted whitelist rule owns -> gene-wide lists skip them
        self.restricted={(r['gene'],c) for r in self.rules if r['action']=='whitelist' and (r['exons'] or r['aa'])
                         for c in r['classes']}

    # checks whether a variant matches a specific, custom rule
    def _rule_hits(self,rule,classes,row):
        if not rule['classes']&classes:
            return False

        if rule['aa']:
            # check if the variant's residue span overlaps a range
            span=_aa_span(row)
            if span is None or not any(lo<=span[1] and span[0]<=hi for lo,hi in rule['aa']):
                return False

        if rule['exons']:
            # if it has an INTRON value
            if 'splice' in classes and row['INTRON'] not in ('','.'):
                n=_region_num(row['INTRON'])          # intron n = between exons n and n+1
                # check if the intron is adjacent to a whitelisted exon
                ok=n is not None and any(lo-1<=n<=hi for lo,hi in rule['exons'])
            else:
                # check if the exon is in the whitelist
                n=_region_num(row['EXON'])
                ok=n is not None and any(lo<=n<=hi for lo,hi in rule['exons'])
            if not ok:
                return False
            
        return True

    def annotate(self,row):
        # initialize whitelist status and matching rules
        out={'wl_whitelist':False,'wl_mis':False,'wl_lof':False,'wl_splice':False,'wl_exception':False,
             'wl_manualreview':False,'wl_rule':[],'wl_transcript':'.'}

        gene=row['Gene']

        if gene in self.accession:
            classes=var_classes(row['Variant.Consequence'])
            mane=(row['MANE.Select'] or '').split('.')[0]

            # TRUE / FALSE = MANE transcript (version dropped) matches the whitelist Accession; '.' = no MANE / Accession
            if mane and self.accession[gene]:
                out['wl_transcript']=mane==self.accession[gene]
            
            # check if missense and protein change is in the missense whitelist
            pc=protein_change(row['Amino.acids'],row['Protein.position'])
            if pc and (gene,pc) in self.missense:
                out['wl_mis']=True
                out['wl_rule'].append('missense')
            
            # check if variant is in a LoF gene
            lof=classes&CLASS_ALIASES['lof']
            if lof and gene in self.lof_genes and not any((gene,c) in self.restricted for c in lof):
                out['wl_lof']=True
                out['wl_rule'].append('lof')

            # check if variant is in a splice gene
            if 'splice' in classes and gene in self.splice_genes and (gene,'splice') not in self.restricted:
                out['wl_splice']=True
                out['wl_rule'].append('splice')
            
            # check if variant matches a custom rule
            for rule in self.rules:
                # go through each rule & skip if it does not match
                if rule['gene']!=gene or not self._rule_hits(rule,classes,row):
                    continue
                
                if rule['action']=='manualreview':
                    out['wl_manualreview']=True
                    out['wl_rule'].append('manual_review')
                    continue
                
                out['wl_exception']=True

                out['wl_rule'].append(f"{rule['label']}")

                hit=rule['classes']&classes

                out['wl_mis']|='missense' in hit
                out['wl_lof']|=bool(hit&CLASS_ALIASES['lof'])
                out['wl_splice']|='splice' in hit
            
            # overall whitelist flag
            out['wl_whitelist']=out['wl_mis'] or out['wl_lof'] or out['wl_splice'] or out['wl_exception']
        
        # format output values
        out['wl_rule']='; '.join(out['wl_rule']) or '.'

        return {k:('TRUE' if v is True else 'FALSE' if v is False else v) for k,v in out.items()}


# ---- tumour hotspots
NEW_COLS=['wl_hotspot','wl_hotspot_rule','protected']


def load_hotspots(path):
    """hotspot_rules.tsv -> rule dicts for Whitelist._rule_hits"""
    rules=[]
    for i,r in enumerate(_rows(path),1):
        classes=set()
        for c in (x.strip() for x in r['var_class'].split(',') if x.strip()):
            if c not in VALID_CLASSES:
                raise SystemExit(f'{path} row {i}: unknown var_class {c!r}')
            classes|=CLASS_ALIASES.get(c,{c})
        if not classes:
            raise SystemExit(f'{path} row {i}: var_class is empty')
        # empty aa_ranges = anywhere in the gene
        label=f"{r['gene']}_{r['var_class']}"+(f"_{r['aa_ranges']}" if r['aa_ranges'] else '')
        rules.append(dict(gene=r['gene'],classes=classes,exons=[],aa=_ranges(r['aa_ranges']),label=label))
    return rules


def main(argv=None):
    p=argparse.ArgumentParser()
    p.add_argument('-i','--input',required=True,help='merged.csv')
    p.add_argument('-o','--output',required=True)
    p.add_argument('--protected',required=True,help='the protected rows only')
    p.add_argument('--stats',required=True)
    p.add_argument('--whitelist-dir',required=True,help='CHIP whitelist files')
    p.add_argument('--hotspots',required=True,help='hotspot_rules.tsv')
    a=p.parse_args(argv)
    wl=Whitelist(a.whitelist_dir)
    rules=load_hotspots(a.hotspots)

    n=Counter()
    with open(a.input,newline='') as fi, open(a.output,'w',newline='') as fo, open(a.protected,'w',newline='') as fp:
        r=csv.DictReader(fi)
        missing=[c for c in REQUIRED if c not in (r.fieldnames or [])]
        if missing:
            raise SystemExit(f'{a.input} is missing columns: {missing}')
        w=csv.DictWriter(fo,fieldnames=r.fieldnames+WL_COLS+NEW_COLS)
        w_pr=csv.DictWriter(fp,fieldnames=r.fieldnames+WL_COLS+NEW_COLS)
        w.writeheader()
        w_pr.writeheader()
        for x in r:
            x.update(wl.annotate(x))
            classes=var_classes(x['Variant.Consequence'])
            hits=[h['label'] for h in rules if h['gene']==x['Gene'] and wl._rule_hits(h,classes,x)]
            x['wl_hotspot']=str(bool(hits)).upper()
            x['wl_hotspot_rule']='; '.join(hits) or '.'
            x['protected']=str('TRUE' in (x['wl_whitelist'],x['wl_manualreview'],x['wl_hotspot'])).upper()
            n['calls']+=1
            for c in ('wl_whitelist','wl_manualreview','wl_hotspot','protected'):
                n[c]+=x[c]=='TRUE'
            w.writerow(x)
            if x['protected']=='TRUE':
                w_pr.writerow(x)

    order=['calls','wl_whitelist','wl_manualreview','wl_hotspot','protected']
    with open(a.stats,'w') as fs:
        fs.write('step\tn\n'+''.join(f'{k}\t{n[k]}\n' for k in order))
    print(', '.join(f'{k}={n[k]}' for k in order)+f' -> {a.output}')
    print('whitelist files: '+', '.join(f'{k}={v}' for k,v in wl.paths.items())+f', hotspots={a.hotspots}')


if __name__=='__main__':
    main()
