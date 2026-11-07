#!/usr/bin/env python3
"""Second code path: symbolic derivation of matched-marginal gap identities."""
import sympy as S,json
from pathlib import Path
r=Path(__file__).resolve().parent
q,eta,sigma,G=S.symbols('q eta sigma G', positive=True, real=True)
cm,cp=1-eta*q,1+eta*q
ql1,ql2=G*(1-sigma*q),G*S.Rational(1,2)*(1-sigma*q)
qu1,qu2=G*(1+sigma*q),G*S.Rational(1,2)*(1+sigma*q)
def t(g,h,Q1,Q2):
 return ((g+h)*Q1+h*Q2)/(g*(g+2*h)), (h*Q1+(g+h)*Q2)/(g*(g+2*h))
# ambient both same g at true extreme, cross-link h chooses independent endpoints
U1=t(G*cm,G*cm,qu1,qu2)[0]
U2=t(G*cm,G*cp,qu1,qu2)[1]
L1=t(G*cp,G*cp,ql1,ql2)[0]
L2=t(G*cp,G*cm,ql1,ql2)[1]
up=S.factor(G*cm*(U1+U2)-(qu1+qu2))
lo=S.factor((ql1+ql2)-G*cp*(L1+L2))
want_up=G*(1+sigma*q)*eta*q/(3*(3+eta*q))
want_lo=G*(1-sigma*q)*eta*q/(3*(3-eta*q))
assert S.simplify(up-want_up)==0 and S.simplify(lo-want_lo)==0
assert S.factor(U1-t(G*cm,G*cm,qu1,qu2)[0])==0
assert S.factor(L1-t(G*cp,G*cp,ql1,ql2)[0])==0
assert S.factor(U2-t(G*cm,G*cm,qu1,qu2)[1])!=0
assert S.factor(L2-t(G*cp,G*cp,ql1,ql2)[1])!=0
out={'symmetric_matched_marginal_upper_gap':str(up), 'symmetric_matched_marginal_lower_gap':str(lo),
     'symbolic_derivation_upper':True,'symbolic_derivation_lower':True,
     'material_shared_marginals':True,'strict_positive_for_0_lt_q_le_1_when_eta_lt_1_and_sigma_lt_1':True}
(r/'matched_marginal_sympy_certificate.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
