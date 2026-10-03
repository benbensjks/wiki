"""Verify the four errata claims against dense_all.csv and perbit.csv."""
import csv, glob, os, sys
sys.stdout.reconfigure(encoding='utf-8')
HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, 'results')


def load(pat, name):
    d = sorted(glob.glob(os.path.join(R, pat)))[-1]
    with open(os.path.join(d, name), encoding='utf-8') as fh:
        return d, list(csv.DictReader(fh))


dd, dense = load('dense_*', 'dense_all.csv')
pd_, perbit = load('perbit_*', 'perbit.csv')
print('dense =', os.path.basename(dd))
print('perbit =', os.path.basename(pd_))
print()

def g(r, k):
    v = r.get(k)
    return None if v in (None, '', 'None') else v

print('=== claim 1+2: upper uM boundary at n=6 (dense) ===')
hdr = ('%-16s %-4s %-4s %-6s %-10s %-10s %-10s' %
       ('label', 'cnt', 'ev', 'full', 'commit', 's0_rev/g/f', 's1_rev/g/f'))
print(hdr); print('-' * len(hdr))
for r in dense:
    if r['family'] != 'DENSE_HIGH' or r['n_A1_gate'] != '6.0':
        continue
    print('%-16s %-4s %-4s %-6s %-10s %-10s %-10s' % (
        r['label'], r['counting_passed'], r['event_causality_passed'],
        r['full_certified'], g(r, 'commitment'),
        '%s/%s/%s' % (g(r, 's0_rev'), g(r, 's0_gate'), g(r, 's0_flip')),
        '%s/%s/%s' % (g(r, 's1_rev'), g(r, 's1_gate'), g(r, 's1_flip'))))
print(' S2 margin_S2_h by uM:')
for r in dense:
    if r['family'] == 'DENSE_HIGH' and r['n_A1_gate'] == '6.0':
        print('    uM=%-6s S2m=%-10s bandM=%s' %
              (r['value'], g(r, 'margin_S2_h'), g(r, 'band_margin_h')))
print()

print('=== claim 3: lower uM boundary, bit0 vs bit1 (dense) ===')
hdr = ('%-16s %-4s %-4s %-10s %-10s %-12s' %
       ('label', 'cnt', 'ev', 's0 r/g/f', 's1 r/g/f', 'sequence'))
print(hdr); print('-' * len(hdr))
for r in dense:
    if r['family'] != 'DENSE_LOW':
        continue
    if r['n_A1_gate'] not in ('6.0',) or r['value'] not in ('4.025', '4.3125', '4.6', '3.7375'):
        continue
    print('%-16s %-4s %-4s %-10s %-10s %-12s' % (
        r['label'], r['counting_passed'], r['event_causality_passed'],
        '%s/%s/%s' % (g(r, 's0_rev'), g(r, 's0_gate'), g(r, 's0_flip')),
        '%s/%s/%s' % (g(r, 's1_rev'), g(r, 's1_gate'), g(r, 's1_flip')),
        r['sequence'][:20]))
print()

print('=== claim 4: K_A[1] boundary data ===')
for fam in ('DENSE_KA', 'DENSE_KF'):
    print(' --', fam)
    for r in dense:
        if r['family'] != fam:
            continue
        print('    n=%-4s value=%-8s f=%-6s cnt=%s full=%s' % (
            r['n_A1_gate'], r['value'], r['value'],
            r['counting_passed'], r['full_certified']))
print()

print('=== n=8 at uM x1.35 (event arm differs?) ===')
for r in dense:
    if r['family'] == 'DENSE_HIGH' and r['value'] == '7.7625':
        print('    n=%-4s cnt=%s ev=%s s1=%s/%s/%s' % (
            r['n_A1_gate'], r['counting_passed'], r['event_causality_passed'],
            g(r, 's1_rev'), g(r, 's1_gate'), g(r, 's1_flip')))
print()

print('=== perbit: S0/S1/S2 by n ===')
for r in perbit:
    if r.get('label', '').startswith('n='):
        print('  %-16s S0: bm=%-9s hole=%-8s | S1: bm=%-9s hold=%-8s | '
              'S2: bm=%-9s hold=%-8s' % (
                  r['label'], g(r, 'S0_band_margin'), g(r, 'S0_hold_h'),
                  g(r, 'S1_band_margin'), g(r, 'S1_hold_h'),
                  g(r, 'S2_band_margin'), g(r, 'S2_hold_h')))
