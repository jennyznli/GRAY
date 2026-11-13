"""
Split a cohort VEP csv (annotate.smk, parser_mode cohort) into the calls / sites / sample_qc inputs of the CHIP scripts.

    python scripts/chip_prep.py -i cohort.vep.csv -o OUTDIR [--samples-from unfiltered.vep.csv]

-i: cohort vep csv, optionally pre-filtered (row_id / n_carriers / sample_freq cols are dropped).
--samples-from: unfiltered csv for the sample list, so samples w/ no remaining calls still count.
"""
import csv
import gzip
import argparse
from pathlib import Path

DROP_COLS = {'row_id', 'n_carriers', 'sample_freq'}  # basic_filter.py extras
SAMPLE_COLS = ['Sample.Zyg', 'Sample.Depth', 'Sample.AltDepth', 'Sample.AltFrac', 'Sample.AF',
               'Sample.GQ', 'Sample.FAD', 'Sample.F1R2', 'Sample.F2R1']


def main():
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--input', required=True)
    p.add_argument('-o', '--outdir', required=True)
    p.add_argument('--samples-from', default='')
    a = p.parse_args()
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)

    sites, calls, samples = {}, {}, set()
    with open(a.input, newline='') as fh:
        r = csv.DictReader(fh)
        site_cols = [c for c in r.fieldnames if c not in SAMPLE_COLS and c not in DROP_COLS and c != 'Sample.ID']
        for x in r:
            s, vid = x['Sample.ID'], x['ID']
            samples.add(s)
            sites.setdefault((vid, x['Feature.Accession']), {c: x[c] for c in site_cols} | {'FILTER': '.'})
            if (s, vid) in calls:  # one call per sample x variant; annotation rows repeat it
                continue
            calls[(s, vid)] = {'Sample.ID': s, 'variant_id': vid, 'FILTER': x['FILTER']} | \
                {c: x.get(c, '') for c in SAMPLE_COLS}

    if a.samples_from:
        with open(a.samples_from, newline='') as fh:
            samples |= {x['Sample.ID'] for x in csv.DictReader(fh)}

    with gzip.open(out / 'calls.csv.gz', 'wt', newline='') as fo:
        w = csv.DictWriter(fo, fieldnames=['Sample.ID', 'variant_id', 'FILTER'] + SAMPLE_COLS)
        w.writeheader()
        w.writerows(calls.values())
    with open(out / 'sites.csv', 'w', newline='') as fo:
        w = csv.DictWriter(fo, fieldnames=site_cols)
        w.writeheader()
        w.writerows(sites.values())
    with open(out / 'sample_qc.csv', 'w', newline='') as fo:
        fo.write('sample_id,sample_flags\n' + ''.join(f'{s},PASS\n' for s in sorted(samples)))
    print(f'{len(samples)} samples, {len(calls)} calls, {len(sites)} site rows -> {out}')


if __name__ == '__main__':
    main()
