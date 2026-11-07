#!/usr/bin/env python3
"""Independent exact Fraction corner enumeration for the Table-I1 descriptive eta scenario.
3232 mathematical grid cases, NOT measured thermal states or a continuum proof.
Continuum proof inherited unchanged from R4; same sign conditions hold since eta<0.01.
"""
from fractions import Fraction as F
from itertools import product
import json, csv
from pathlib import Path
P=Path(__file__).resolve().parent
src=json.loads((P/'NIST15_source_audit_certificate.json').read_text())
e=F(src['max_abs_deviation_exact_fraction']);assert 0<e<F(1,100)
g0=F(331467,1016000); s=F(1,5)

def solve(g1,g2,h,Q1,Q2):
    det=g1*g2+h*(g1+g2)
    assert det>0
    t1=((g2+h)*Q1+h*Q2)/det
    t2=(h*Q1+(g1+h)*Q2)/det
    assert g1*t1+h*(t1-t2)==Q1
    assert g2*t2+h*(t2-t1)==Q2
    assert g1*t1+g2*t2==Q1+Q2
    return t1,t2

summary=[];cases=0
for j in range(101):
    q=F(j,100)
    lo=g0*(1-e*q);hi=g0*(1+e*q)
    q1lo=g0*(1-s*q);q1hi=g0*(1+s*q)
    q2lo=g0*F(1,2)*(1-s*q);q2hi=g0*F(1,2)*(1+s*q)
    temps=[]
    for g1,g2,h in product((lo,hi),repeat=3):
      for Q1,Q2 in product(((q1lo,q1hi)),((q2lo,q2hi))):
        t=solve(g1,g2,h,Q1,Q2);temps.append(t);cases+=1
    minT=[min(t[i] for t in temps) for i in (0,1)]
    maxT=[max(t[i] for t in temps) for i in (0,1)]
    U1_expected=solve(lo,lo,lo,q1hi,q2hi)[0]
    U2_expected=solve(lo,lo,hi,q1hi,q2hi)[1]
    L1_expected=solve(hi,hi,hi,q1lo,q2lo)[0]
    L2_expected=solve(hi,hi,lo,q1lo,q2lo)[1]
    assert maxT == [U1_expected,U2_expected]
    assert minT == [L1_expected,L2_expected]
    upper_extra=lo*sum(maxT)-q1hi-q2hi
    lower_deficit=q1lo+q2lo-hi*sum(minT)
    expected_upper=g0*(1+s*q)*e*q/(3*(3+e*q))
    expected_lower=g0*(1-s*q)*e*q/(3*(3-e*q))
    assert upper_extra == expected_upper
    assert lower_deficit == expected_lower
    assert upper_extra >= 0 and lower_deficit>=0
    # Shared-multiplier extreme control must be attained by the SAME complete vector
    shared_hi=solve(lo,lo,lo,q1hi,q2hi)
    shared_lo=solve(hi,hi,hi,q1lo,q2lo)
    assert shared_hi[0] == maxT[0]
    assert shared_lo[0] == minT[0]
    summary.append({'alpha':str(1-q),'q':str(q),
     'upper_excess_W':f'{float(upper_extra):.12g}',
     'lower_deficit_W':f'{float(lower_deficit):.12g}',
     'upper_L1_potential_K':f'{float(upper_extra/lo):.12g}',
     'lower_L1_potential_K':f'{float(lower_deficit/hi):.12g}'})
assert cases==3232
# monotonic nesting: alpha increases (q decreases) -> lower nondecreasing, upper nonincreasing
for q in [F(k,100) for k in range(101)]:
 assert g0*(1-e*q)>0
# independent derivative signs established already in R4 for e<=1/100 and sigma=1/5.
with (P/'NIST15_descriptive_eta_101cut_grid.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(summary[0]));w.writeheader();w.writerows(summary)
result={'source_inputs':'NIST Table I1, DOI 10.6028/jres.119.012; certified fit DOI 10.6028/NIST.SP.260-173',
        'max_observed_residual_fraction':src['max_abs_deviation_exact_fraction'],
        'grid_levels':101,'corner_sources_per_level':32,'exact_conservative_corner_states':cases,
        'observed_support_upper_energy_gap_W':summary[-1]['upper_excess_W'],
        'observed_support_lower_energy_gap_W':summary[-1]['lower_deficit_W'],
        'proof_status':'The R4 all-alpha analytic proof applies under declared eta<0.01; grid tests are diagnostics only',
        'empirical_status':'historical published material measurements only, not physical-network validation'}
(P/'NIST15_independent_corner_certificate.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
