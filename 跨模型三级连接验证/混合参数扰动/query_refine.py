import glob, json, os, sys
sys.stdout.reconfigure(encoding='utf-8')
d = sorted(glob.glob(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                  'results', 'refine_*')))[-1]
j = json.load(open(os.path.join(d, 'refine_all.json'), encoding='utf-8'))
print('dir', os.path.basename(d))
print()
hdr = ('%-30s %-6s %-5s %-6s %-6s %-6s %-6s %-6s %-6s' %
       ('label', 'cert', 'S2m', 's0', 's1', 'o2o', 'ord', 'alt', 'inc'))
print(hdr)
print('-' * len(hdr))
for r in j['rows']:
    s0 = r.get('stage0') or {}
    s1 = r.get('stage1') or {}
    a = s0.get('one_to_one'); b = s0.get('order'); c = s0.get('alt')
    e = s1.get('one_to_one'); f = s1.get('order'); g = s1.get('alt')
    m = r.get('margin_S2_h')
    print('%-30s %-6s %-5s %-6s %-6s %-6s %-6s %-6s %-6s' % (
        r['label'][:30], r['certified'],
        'None' if m is None else '%.3f' % m,
        '%d/%d/%d' % (s0.get('reverse') or 0, s0.get('gate') or 0, s0.get('flip') or 0),
        '%d/%d/%d' % (s1.get('reverse') or 0, s1.get('gate') or 0, s1.get('flip') or 0),
        'T' if a else 'F', 'T' if b else 'F', 'T' if c else 'F',
        'T' if r['increments_mod8'] else 'F'))
print()
print('--- failure-mode detail (non-certified rows) ---')
for r in j['rows']:
    if r['certified']:
        continue
    s0 = r.get('stage0') or {}
    s1 = r.get('stage1') or {}
    print('*', r['label'])
    print('   seq      :', r['sequence'])
    print('   inc_mod8 :', r['increments_mod8'], '| clips', r['boundary_clips'],
          '| min_commit', r['minimum_commitment'])
    print('   stage0   :', {k: s0.get(k) for k in
                            ('reverse', 'gate', 'flip', 'one_to_one', 'order', 'alt')})
    print('   stage1   :', {k: s1.get(k) for k in
                            ('reverse', 'gate', 'flip', 'one_to_one', 'order', 'alt')})
    print('   clock    :', (r.get('signal_ranges') or {}).get('clock'))
