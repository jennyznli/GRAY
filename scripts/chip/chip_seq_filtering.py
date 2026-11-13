"""
CHIP sequence-level filtering for one chromosome: join somatic_variants calls to their VEP annotation and apply CH specific filters. 

Steps:
  01  Gene in Full_CHIP_gene_list_*.txt, Bio.type == protein_coding,
      Sample.AltFrac / Sample.Depth / Sample.AltDepth all present (not '', '.', 'NA', '-1')
  02  MANE.Select or MANE.PlusClinical filled
  03  SOMQC has no 'orientation' flag (filter.smk: alt F1R2 and F2R1 >= 1); --no-orientation skips it
  04  Sample.Depth >= --min-dp
  05  Sample.AltDepth >= --min-alt-depth
  06  Sample.AltFrac > --min-vaf (strict)
Mutect2 FILTER is carried through, not used.

Output: one row per sample x variant x CH gene. Columns = Sample.ID, the annotate csv's columns (its
'.' FILTER replaced by the sample's FILTER), then SOMQC + the Sample.* columns. Default = rows passing all
steps; --keep-all writes every CH-gene row with chip_s01..chip_s06 TRUE/FALSE columns.

--stats: long tsv  chrom(--chrom)  step  n   (input = all calls on the chromosome; then rows passing each step)
"""
import os
import sys
import csv
import gzip
import argparse
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chip_whitelist import _latest, _rows  # noqa: E402

STEPS = ['s01_gene_biotype_values', 's02_mane', 's03_orientation', 's04_dp', 's05_alt_depth', 's06_vaf']
STEP_COLS = [f'chip_{s[:3]}' for s in STEPS]
SAMPLE_COLS = ['SOMQC', 'Sample.Depth', 'Sample.AltDepth', 'Sample.AltFrac', 'Sample.AF',
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
    p.add_argument('--calls', required=True, help='somatic_variants.{chrom}.csv.gz')
    p.add_argument('--sites', required=True, help='annotate csv for the same chrom')
    p.add_argument('-o', '--output', required=True)
    p.add_argument('--stats', required=True)
    p.add_argument('--chrom', required=True, help='label for the --stats rows')
    p.add_argument('--whitelist-dir', default='whitelist')
    p.add_argument('--min-dp', type=int, default=20)
    p.add_argument('--min-alt-depth', type=int, default=3)
    p.add_argument('--min-vaf', type=float, default=0.02)
    p.add_argument('--no-orientation', action='store_true')
    p.add_argument('--keep-all', action='store_true')
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
    out_cols += SAMPLE_COLS + (STEP_COLS if a.keep_all else [])

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
                      a.no_orientation or 'orientation' not in c['SOMQC'].split(';'),
                      dp is not None and dp >= a.min_dp,
                      ad is not None and ad >= a.min_alt_depth,
                      vaf is not None and vaf > a.min_vaf]
                n['in_gene'] += 1
                cum = True
                for step, s_ok in zip(STEPS, ok):
                    cum = cum and s_ok
                    n[step] += cum
                if cum or a.keep_all:
                    row = site | {'Sample.ID': c['Sample.ID'], 'FILTER': c['FILTER']}
                    row.update({k: c[k] for k in SAMPLE_COLS})
                    if a.keep_all:
                        row.update({col: 'TRUE' if s_ok else 'FALSE' for col, s_ok in zip(STEP_COLS, ok)})
                    w.writerow(row)

    order = ['input', 'in_gene'] + STEPS
    with open(a.stats, 'w') as fs:
        fs.write('chrom\tstep\tn\n')
        for k in order:
            fs.write(f'{a.chrom}\t{k}\t{n[k]}\n')
    print(f'{a.calls}: ' + ', '.join(f'{k}={n[k]}' for k in order) + f' -> {a.output}')
    print(f'min_dp={a.min_dp} min_alt_depth={a.min_alt_depth} min_vaf>{a.min_vaf} '
          f'orientation={"skipped" if a.no_orientation else "required"}; {len(genes)} genes in {gene_file}')


if __name__ == '__main__':
    main()
