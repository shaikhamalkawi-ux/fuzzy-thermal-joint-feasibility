#!/usr/bin/env python3
"""Exact Fraction certificate for physical-scale material-anchored COMMON conductivity.
   This is a synthetic geometry/loads example, not an experimental validation.
   NIST source: https://doi.org/10.6028/NIST.SP.260-173
"""
import csv, json, hashlib, itertools, os
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parent
K_REF=F(110489,10**9)*300 # 1.10489e-4 * 300 W/(m K)
A=F(1,4) # 0.25 m^2 chosen scenario, not NIST specimen observation
L=F(127,5000) # 0.0254 m scenario using NIST nominal 25.4mm format
C=A/L
G_REF=K_REF*C
ETA=F(1,100) # 1% illustrative fuzzy-scale envelope; NOT NIST statistical confidence interval
SIGMA=F(1,5) # 20% hypothetical input uncertainty; not NIST measurement
assert K_REF == F(331467,10**7)

def solve(M,b):
    n=len(M)
    a=[list(M[i])+[b[i]] for i in range(n)]
    for j in range(n):
        p=next((i for i in range(j,n) if a[i][j]), None)
        assert p is not None
        a[j],a[p]=a[p],a[j]
        d=a[j][j];a[j]=[v/d for v in a[j]]
        for i in range(n):
            if i!=j:
                f=a[i][j]; a[i]=[a[i][k]-f*a[j][k] for k in range(n+1)]
    return [a[i][-1] for i in range(n)]

def mv(M,x):return [sum(row[j]*x[j] for j in range(len(x))) for row in M]

def inv(H):
    n=len(H)
    return [[solve(H,[F(int(i==j)) for i in range(n)])[r] for j in range(n)] for r in range(n)]

def fmt(x):return str(x.numerator) if x.denominator==1 else str(x.numerator)+"/"+str(x.denominator)

CASE={
 "2node": {"H":[[F(2), F(-1)],[F(-1), F(2)]],"v":[F(1),F(1,2)]},
 "3node": {"H":[[F(2), F(-1),F(0)],[F(-1),F(3),F(-1)],[F(0),F(-1),F(2)]],"v":[F(1),F(1,2),F(1,5)]}
}
summary={"source":"NIST SP 260-173 10.6028/NIST.SP.260-173", "source_k_at_300_W_mK":fmt(K_REF),
 "source_expanded_U_percent":1, "scenario_geometry_A_m2":fmt(A),"scenario_geometry_L_m":fmt(L),
 "G_reference_W_K":fmt(G_REF),"scenario_conductivity_support_fraction":fmt(ETA),
 "scenario_heat_support_fraction":fmt(SIGMA),"NIST_1pct_to_fuzzy":"NO: numerical scale only, separate modeling assumption",
 "cases":{},"corner_checks":0,"global_energy_checks":0,"component_extrema_checks":0,"nesting_checks":0}
rows=[]
for case,scenario in CASE.items():
    H=scenario["H"];v=scenario["v"];n=len(v)
    B=inv(H);assert all(a>0 for row in B for a in row),"strictly positive M-matrix inverse required"
    nominal=mv(B,v)
    prev=None
    for k in range(101):
        q=F(k,100);alpha=1-q
        clow=1-ETA*q;chigh=1+ETA*q
        qlow=[G_REF*t*(1-SIGMA*q) for t in v]
        qhigh=[G_REF*t*(1+SIGMA*q) for t in v]
        low=[nominal[i]*(1-SIGMA*q)/chigh for i in range(n)]
        high=[nominal[i]*(1+SIGMA*q)/clow for i in range(n)]
        actual_min=[None]*n;actual_max=[None]*n
        for c in [clow,chigh]:
            K=[[G_REF*c*H[i][j] for j in range(n)] for i in range(n)]
            for selectors in itertools.product((0,1),repeat=n):
                qvector=[qlow[i] if selectors[i]==0 else qhigh[i] for i in range(n)]
                t=solve(K,qvector)
                assert mv(K,t)==qvector
                summary["corner_checks"]+=1
                assert sum(t)*G_REF*c==sum(qvector),"global energy closure"
                summary["global_energy_checks"]+=1
                for i in range(n):
                    actual_min[i]=t[i] if actual_min[i] is None or t[i]<actual_min[i] else actual_min[i]
                    actual_max[i]=t[i] if actual_max[i] is None or t[i]>actual_max[i] else actual_max[i]
        assert actual_min==low and actual_max==high
        summary["component_extrema_checks"]+=2*n
        for c,Q,T in [(clow,qhigh,high),(chigh,qlow,low)]:
            K=[[G_REF*c*H[i][j] for j in range(n)] for i in range(n)]
            assert mv(K,T)==Q,"same network parameter and source attains simultaneous bounding vector"
            assert sum(T)*G_REF*c==sum(Q)
        if prev is not None:
            lastlow,lasthigh=prev
            assert all(low[i]<=lastlow[i] and high[i]>=lasthigh[i] for i in range(n)),"nested cuts"
            summary["nesting_checks"]+=1
        prev=(low,high)
        rows.append([case,fmt(alpha),fmt(q),fmt(clow),fmt(chigh),
           ';'.join(map(fmt,low)),';'.join(map(fmt,high)),
           fmt(sum(qlow)),fmt(sum(qhigh)),"0","0","PASS"])
    q=F(1)
    lo=[nominal[i]*(1-SIGMA)/(1+ETA) for i in range(n)]
    hi=[nominal[i]*(1+SIGMA)/(1-ETA) for i in range(n)]
    summary["cases"][case]={"nominal_dT_K":list(map(fmt,nominal)),"support_dT_low_K":list(map(fmt,lo)),
       "support_dT_high_K":list(map(fmt,hi)),"support_T_absolute_range_K":[[fmt(F(300)+lo[i]),fmt(F(300)+hi[i])] for i in range(n)],
       "independent_input_corners":2**(n+1),"box_high_low_energy_gaps_W":"0 with the appropriate common-k realization"}

with (ROOT/'NIST_commonK_101cut_certificate.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.writer(f);w.writerow(['case','alpha','q','conductivity_multiplier_low','conductivity_multiplier_high',
       'temperature_lower_dT_K','temperature_upper_dT_K','Q_sum_min_W','Q_sum_max_W',
       'upper_joint_energy_gap_W','lower_joint_energy_gap_W','check']);w.writerows(rows)
with (ROOT/'NIST_commonK_certificate.json').open('w',encoding='utf-8') as f:json.dump(summary,f,indent=2,sort_keys=True)
assert summary['corner_checks']==(8+16)*101==2424
assert summary['nesting_checks']==200
print(json.dumps(summary,indent=2))
