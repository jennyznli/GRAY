"""
Split a gray cohort VEP csv into the calls / sites / sample_qc inputs of the CHIP scripts (copied from 00-pmbb).

    python scripts/chip/chip_example_prep.py -i all_samples.mutect2.no_pon.vep.csv -o OUTDIR [--min-dp 20 --min-ad 3]

OUTDIR/calls.csv.gz   Sample.ID, variant_id, FILTER, SOMQC (somatic_flags.py rules, alt-only orientation), Sample.*
OUTDIR/sites.csv      one annotation row per ID x Feature.Accession (FILTER = '.')
OUTDIR/sample_qc.csv  every sample PASS (gray has no qc.smk metrics)
"""
import csv
import gzip
import argparse
from pathlib import Path

SAMPLE_COLS = ['Sample.Zyg', 'Sample.Depth', 'Sample.AltDepth', 'Sample.AltFrac', 'Sample.AF',
               'Sample.GQ', 'Sample.FAD', 'Sample.F1R2', 'Sample.F2R1']


def alt(x):
    try:
        return int(x.split('|')[1])
    except (AttributeError, IndexError, ValueError):
        return None


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def main():
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--input', required=True)
    p.add_argument('-o', '--outdir', required=True)
    p.add_argument('--min-dp', type=int, default=20)
    p.add_argument('--min-ad', type=int, default=3)
    a = p.parse_args()
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)

    sites, calls, samples = {}, {}, set()
    with open(a.input, newline='') as fh:
        r = csv.DictReader(fh)
        site_cols = [c for c in r.fieldnames if c not in SAMPLE_COLS and c != 'Sample.ID']
        for x in r:
            s, vid = x['Sample.ID'], x['ID']
            samples.add(s)
            sites.setdefault((vid, x['Feature.Accession']), {c: x[c] for c in site_cols} | {'FILTER': '.'})
            if (s, vid) in calls:  # one call per sample x variant; annotation rows repeat it
                continue
            dp, ad = num(x['Sample.Depth']), num(x['Sample.AltDepth'])
            o1, o2 = alt(x.get('Sample.F1R2')), alt(x.get('Sample.F2R1'))  # VarDict: absent -> orientation flag
            flags = [f for f, bad in [('low_dp', dp is None or dp < a.min_dp),
                                      ('low_ad', ad is None or ad < a.min_ad),
                                      ('orientation', o1 is None or o2 is None or o1 < 1 or o2 < 1)] if bad]
            calls[(s, vid)] = {'Sample.ID': s, 'variant_id': vid, 'FILTER': x['FILTER'],
                               'SOMQC': ';'.join(flags) or 'PASS'} | {c: x.get(c, '') for c in SAMPLE_COLS}

    with gzip.open(out / 'calls.csv.gz', 'wt', newline='') as fo:
        w = csv.DictWriter(fo, fieldnames=['Sample.ID', 'variant_id', 'FILTER', 'SOMQC'] + SAMPLE_COLS)
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
