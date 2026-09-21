"""debug_conservation.py — locate the DNA conservation leak numerically."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
import zhao_core as zc

P = zc.default_params()
C = zc._rate_constants(P)

def int_src(t):
    return P['k_int'] * zc.square_pulse(t, P)

# random-ish positive state
rng = np.random.default_rng(0)
y = np.abs(rng.normal(0, 0.05, 38))
F = zc.rhs(0.3, y, P, C, int_src)

dna_idx = zc._LR_IDX + zc._PB_IDX
print("sum d(DNA)/dt =", F[dna_idx].sum())

# per-species table of DNA species derivatives
names = {0:'LR',4:'PB-int2',5:'LR-int2',6:'PB-int4',7:'LR-int4',8:'LR-ints2',
         9:'PB-ints',10:'LR-ints1',11:'PB-int2-rdf2',12:'LR-int2-rdf2',
         13:'PB-int4-rdf4',14:'LR-int4-rdf4',15:'PB-int-rdfs2',16:'PB-int-rdfs1',
         17:'LR-int-rdfs',20:'PB-int2-rdf',21:'LR-int2-rdf',23:'PB',
         24:'PB-int4-rdf',25:'PB-int4-rdf2',26:'PB-int4-rdf3',27:'LR-int4-rdf',
         28:'LR-int4-rdf2',29:'LR-int4-rdf3',30:'PB-int6i',31:'PB-int6-rdf4i',
         33:'PB-int6-rdfi',34:'PB-int6-rdf2i',35:'PB-int6-rdf3i'}
for i in sorted(names):
    print(f"  y[{i:2d}] {names[i]:16s} F={F[i]:+.6e}")
