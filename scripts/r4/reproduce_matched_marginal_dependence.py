#!/usr/bin/env python3
"""Exact same-marginal uncertainty control: 2-node material conductance links.
Common factor c versus independent g1,g2,h all in same [1±1%*q] band.
Illustrative NIST-property scale; NO calibration/experiment.
"""
import csv,json,itertools
from fractions import Fraction as F
from pathlib import Path
R=Path(__file__).resolve().parent
G=F(331467,1016000) # NIST conductivity at 300K * chosen A/L
eta=F(1,100);sig=F(1,5)

def f(x):return str(x)

def solve2(g1,g2,h,q1,q2):
 d=g1*g2+h*(g1+g2)
 return ((g2+h)*q1+h*q2)/d,(h*q1+(g1+h)*q2)/d

n=0; same_matched=0; rows=[]
for j in range(101):
 q=F(j,100); cmin=1-eta*q;cmax=1+eta*q
 Qlo=(G*(1-sig*q),G*F(1,2)*(1-sig*q))
 Qhi=(G*(1+sig*q),G*F(1,2)*(1+sig*q))
 states=[]
 for cm in itertools.product([cmin,cmax], repeat=3):
  g1,g2,h=[G*x for x in cm]
  for bb in itertools.product([0,1],repeat=2):
   qi=(Qlo[0] if bb[0]==0 else Qhi[0],Qlo[1] if bb[1]==0 else Qhi[1])
   u=solve2(g1,g2,h,*qi)
   assert g1*u[0]+g2*u[1]==qi[0]+qi[1]
   n+=1
   states.append((cm,bb,u))
 M=[min(s[2][i] for s in states) for i in range(2)]
 U=[max(s[2][i] for s in states) for i in range(2)]
 # Exact global extrema identified by derivative sign / positive inverse
 L1=solve2(G*cmax,G*cmax,G*cmax,*Qlo)[0]
 L2=solve2(G*cmax,G*cmax,G*cmin,*Qlo)[1]
 U1=solve2(G*cmin,G*cmin,G*cmin,*Qhi)[0]
 U2=solve2(G*cmin,G*cmin,G*cmax,*Qhi)[1]
 assert M==[L1,L2] and U==[U1,U2]
 same_matched+=1
 # A necessary total energy condition for any independent material network
 # is G*cmin*sum(U)<=sum(Qhi), G*cmax*sum(L)>=sum(Qlo) if U/L attainable.
 up=G*cmin*(U1+U2)-sum(Qhi)
 lo=sum(Qlo)-G*cmax*(L1+L2)
 formula_up=G*(1+sig*q)*eta*q/(3*(3+eta*q))
 formula_lo=G*(1-sig*q)*eta*q/(3*(3-eta*q))
 assert up==formula_up and lo==formula_lo and up>=0 and lo>=0
 # Common factor extremes match node1, and node2 differs for q>0
 CL=solve2(G*cmax,G*cmax,G*cmax,*Qlo)
 CU=solve2(G*cmin,G*cmin,G*cmin,*Qhi)
 assert L1==CL[0] and U1==CU[0]
 assert L2<=CL[1] and U2>=CU[1]
 assert (up>0)==(q>0) and (lo>0)==(q>0)
 rows.append([str(1-q),str(q),str(L1),str(U1),str(L2),str(U2),
              str(CL[0]),str(CU[0]),str(CL[1]),str(CU[1]),str(up),str(lo)])
assert n==3232 and same_matched==101
with (R/'matched_marginal_101cuts.csv').open('w',newline='') as f:
 w=csv.writer(f);w.writerow(['alpha','q','ind_L1_K','ind_U1_K','ind_L2_K','ind_U2_K',
 'common_L1_K','common_U1_K','common_L2_K','common_U2_K','ind_upper_energy_gap_W','ind_lower_energy_gap_W']);w.writerows(rows)
out={'material_conductance_nominal_W_K':str(G),'eta':str(eta),'load_band':str(sig),
     'mathematical_input_marginal_match':'TRUE (same each-link cuts; different dependence)',
     'certificate':'q in [0,1]: ΔU=G(1+σq)ηq/[3(3+ηq)], ΔL=G(1−σq)ηq/[3(3−ηq)]',
     'support_upper_independent_gap_W':str(F(rows[-1][-2])),
     'support_lower_independent_gap_W':str(F(rows[-1][-1])),
     'support_upper_decimal':float(F(rows[-1][-2])),
     'support_lower_decimal':float(F(rows[-1][-1])),
     'support_common_gap_W':'0',
     'exact_independent_corner_solves':n,'exact_global_extrema_levels':same_matched}
(R/'matched_marginal_certificate.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
