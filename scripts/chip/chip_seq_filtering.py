"""
CHIP sequence-level filtering, join chip_prep.py calls to their VEP annotation and apply CH specific filters. 

Steps:
  01  Gene in Full_CHIP_gene_list_*.txt, Bio.type == protein_coding,
      Sample.AltFrac / Sample.Depth / Sample.AltDepth all present (not '', '.', 'NA', '-1')
  02  MANE.Select or MANE.PlusClinical filled
  03  Sample.Depth >= --min-dp
  04  Sample.AltDepth >= --min-alt-depth
  05  Sample.AltFrac > --min-vaf (strict)

Output: one row per sample x variant x CH gene passing all steps. Columns = Sample.ID, the annotate csv's
columns (its '.' FILTER replaced by the sample's FILTER), then the Sample.* columns.

--stats: tsv  step  n   (input = all calls; then rows passing each step)
"""
import os
import sys
import csv
import gzip
import argparse
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chip_whitelist import _latest, _rows  # noqa: E402

STEPS = ['s01_gene_biotype_values', 's02_mane', 's03_dp', 's04_alt_depth', 's05_vaf']
SAMPLE_COLS = ['Sample.Depth', 'Sample.AltDepth', 'Sample.AltFrac', 'Sample.AF',
               'Sample.F1R2', 'Sample.F2R1']
MISSING = ('', '.', 'NA', '-1')


def to_float(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def filled(x):
    return (x or '').strip() not in ('', '.')


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument('--calls', required=True, help='chip_prep.py calls.csv.gz')
    p.add_argument('--sites', required=True, help='chip_prep.py sites.csv')
    p.add_argument('-o', '--output', required=True)
    p.add_argument('--stats', required=True)
    p.add_argument('--whitelist-dir', default='whitelist')
    p.add_argument('--min-dp', type=int, default=20)
    p.add_argument('--min-alt-depth', type=int, default=3)
    p.add_argument('--min-vaf', type=float, default=0.02)
    a = p.parse_args(argv)

    gene_file = _latest(a.whitelist_dir, 'Full_CHIP_gene_list_*.txt')
    genes = {r['Gene'] for r in _rows(gene_file)}

    # CH-gene annotation rows only: variant_id -> [(row, biotype_ok, mane_ok), ...]
    sites = {}
    with open(a.sites, newline='') as fs:
        reader = csv.DictReader(fs)
        site_cols = [c for c in reader.fieldnames if c != 'FILTER']
        for r in reader:
            if r['Gene'] in genes:
                sites.setdefault(r['ID'], []).append(
                    (r, r['Bio.type'] == 'protein_coding',
                     filled(r['MANE.Select']) or filled(r['MANE.PlusClinical'])))

    out_cols = ['Sample.ID'] + site_cols[:]
    # put the sample's FILTER where the parser had it
    i_id = out_cols.index('ID')
    out_cols.insert(i_id + 1, 'FILTER')
    out_cols += SAMPLE_COLS

    n = Counter()
    with gzip.open(a.calls, 'rt', newline='') as fc, open(a.output, 'w', newline='') as fo:
        w = csv.DictWriter(fo, fieldnames=out_cols, extrasaction='ignore')
        w.writeheader()
        for c in csv.DictReader(fc):
            n['input'] += 1
            for site, biotype_ok, mane_ok in sites.get(c['variant_id'], ()):
                dp, ad, vaf = (to_float(c[k]) for k in ('Sample.Depth', 'Sample.AltDepth', 'Sample.AltFrac'))
                ok = [biotype_ok and all(c[k].strip() not in MISSING for k in
                                         ('Sample.AltFrac', 'Sample.Depth', 'Sample.AltDepth')),
                      mane_ok,
                      dp is not None and dp >= a.min_dp,
                      ad is not None and ad >= a.min_alt_depth,
                      vaf is not None and vaf > a.min_vaf]
                n['in_gene'] += 1
                cum = True
                for step, s_ok in zip(STEPS, ok):
                    cum = cum and s_ok
                    n[step] += cum
                if cum:
                    row = site | {'Sample.ID': c['Sample.ID'], 'FILTER': c['FILTER']}
                    row.update({k: c[k] for k in SAMPLE_COLS})
                    w.writerow(row)

    order = ['input', 'in_gene'] + STEPS
    with open(a.stats, 'w') as fs:
        fs.write('step\tn\n')
        for k in order:
            fs.write(f'{k}\t{n[k]}\n')
    print(f'{a.calls}: ' + ', '.join(f'{k}={n[k]}' for k in order) + f' -> {a.output}')
    print(f'min_dp={a.min_dp} min_alt_depth={a.min_alt_depth} min_vaf>{a.min_vaf}; '
          f'{len(genes)} genes in {gene_file}')


if __name__ == '__main__':
    main()
