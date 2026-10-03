import glob, json, os, sys
sys.stdout.reconfigure(encoding='utf-8')
HERE = os.path.dirname(os.path.abspath(__file__))
d = sorted(glob.glob(os.path.join(HERE, 'results', 'swap_*')))[-1]
j = json.load(open(os.path.join(d, 'swap_all.json'), encoding='utf-8'))
print('dir', os.path.basename(d))
print('control (zmh n=6 vs frozen):', j['control_checks'])
print('confound:', j['confound_note'])
print()
rows = {(r['middle'], r['n_A1_gate']): r for r in j['rows']}
keys = ['PB0_high', 'PB0_low', 'PB1_high', 'PB1_low', 'Int1_peak',
        'Int1_source_peak']
sub = ['A0', 'F0', 'Int1', 'T1', 'RDF1', 'Int1_source']
print('=== interface: peak / median_in / median_out ===')
for k in keys:
    line = '%-20s' % k
    for mid in ('zmh', 'm51'):
        for n in (4.0, 6.0):
            v = rows[(mid, n)]['interface'].get(k)
            line += ' %s n=%g:%-12s' % (mid, n, 'None' if v is None else '%.4g' % v)
    print(line)
print()
for s in sub:
    print('--- %s ---' % s)
    for mid in ('zmh', 'm51'):
        for n in (4.0, 6.0):
            v = rows[(mid, n)]['interface'][s]
            print('   %s n=%g  peak=%-10.4g  in=%-10.4g  out=%-10.4g'
                  % (mid, n, v['peak'], v['median_in'], v['median_out']))
print()
print('=== tail dose / shape ===')
for mid in ('zmh', 'm51'):
    for n in (4.0, 6.0):
        r = rows[(mid, n)]
        g = r['g1_window']
        print('  %s n=%g  full=%d  seq=%s' % (mid, n, r['full_certified'], r['sequence']))
        print('      g1 peak=%.4f width_max=%.4f integral=[%.4g,%.4g] far_dwell=%.4g'
              % (r['g1_peak'], g['width_max'], g['integral_min'],
                 g['integral_max'], g['far_dwell_median']))
        print('      Int2 peak=%.3f  J_fwd2=%.4g  J_rev2=%.4g  J_rev1=%.4g'
              % (r['int2_peak'], r['j_fwd2_integral'], r['j_rev2_integral'],
                 r['j_rev1_integral']))
        print('      S2 margin=%s  clock_peaks=%d  s1 rev/gate/flip=%d/%d/%d'
              % (r['margin_S2_h'], r['clock_peak_count'],
                 r['s1']['reverse_events'], r['s1']['gate_events'],
                 r['s1']['flip_events']))
