"""Long confirmation of the sole late-alternating coarse-grid candidate."""
import json
from dataclasses import asdict

from bit0_diagnostic import ROOT, diagnose_one, test_twobit
from model import Extension

out=ROOT/'bit0_results'/'confirmed_candidate'
out.mkdir(parents=True,exist_ok=True)
cfg=asdict(Extension())
cfg.update(uM_per_au=10.0,maturation_half_life_min=20.0,
           complex_on_au_inv_h=0.1,complex_off_h=1.0,add_growth=False)
_,_,cycles,score=diagnose_one(cfg,out,'bit0_180h',180.0)
two=test_twobit(cfg,out,180.0) if score['stable_alternation'] else None
summary=dict(configuration=cfg,bit0_score=score,cycles=int(len(cycles)),twobit=two)
(out/'confirmation_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(summary,ensure_ascii=False,indent=2))
