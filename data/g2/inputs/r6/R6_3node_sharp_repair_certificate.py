#!/usr/bin/env python3
"""Exact symbolic 3-node L1/Linf repair certificate (99fuzzy ICMPC R6).

This verification uses the unchanged R5 3-node synthetic thermal network.
NO new empirical material measurements or solver-ranking claims.
Requires sympy. Python standard library for CSV, json, Fraction.
"""
from __future__ import annotations
import csv, json
from fractions import Fraction
from pathlib import Path
import sympy as s

OUT = Path(__file__).resolve().parent
q = s.symbols('q', real=True)
a,b = s.symbols('a b', positive=True)
v = [s.Integer(2), s.Integer(1), s.Rational(2,5)]
D = 3*a*b+2*a+2*b+1
F = [
 (17*a*b+15*a+20*b+10)/(5*D),
 (17*a*b+15*a+7*b+5)/(5*D),
 (17*a*b+4*a+7*b+2)/(5*D)
]
hm, hp = 1-q/5, 1+q/5
U = [s.factor((5+q)*f.subs({a:x,b:y})) for f,x,y in zip(F,[hm,hp,hp],[hm,hm,hp])]
L = [s.factor((5-q)*f.subs({a:x,b:y})) for f,x,y in zip(F,[hp,hm,hm],[hp,hp,hm])]
Qp, Qm = [[(5+z*q)*w for w in v] for z in (1,-1)]
Dp = s.factor(sum(U)-sum(Qp))
Dm = s.factor(sum(Qm)-sum(L))
Up=[s.factor(u-Dp/3) for u in U]
Lp=[s.factor(l+Dm/3) for l in L]

# Exact inverse reconstruction via local fluxes, not numerical optimization.
def reconstructed_links(t, source):
    ah = s.factor((source[0]-t[0])/(t[0]-t[1]))
    bh = s.factor((t[2]-source[2])/(t[1]-t[2]))
    return ah,bh

apstar,bpstar = reconstructed_links(Up,Qp)
amstar,bmstar = reconstructed_links(Lp,Qm)

# Conservation and local flux equations for each candidate, exact identities.
for label,t,src,aa,bb in [('upper',Up,Qp,apstar,bpstar),('lower',Lp,Qm,amstar,bmstar)]:
    K=s.Matrix([[1+aa,-aa,0],[-aa,1+aa+bb,-bb],[0,-bb,1+bb]])
    assert all(s.cancel(ri)==0 for ri in K*s.Matrix(t)-s.Matrix(src)), label
    assert s.factor(sum(t)-sum(src)) == 0, label

# Prove each reconstructed link is in [1-q/5, 1+q/5] for all q in [0,1].
# Difference sign certificate: numerator = q*P(q), denominator of invariant sign;
# nonnegative coefficients lower-bound for polynomial on [0,1]:
# P(q) >= P(0) + sum_{k>=1} min(0,P_k) > 0, since 0<=q^k<=1.
ineqs = {
 'upper_a_minus_lower':apstar-hm, 'upper_upper_minus_a':hp-apstar,
 'upper_b_minus_lower':bpstar-hm, 'upper_upper_minus_b':hp-bpstar,
 'lower_a_minus_lower':amstar-hm, 'lower_upper_minus_a':hp-amstar,
 'lower_b_minus_lower':bmstar-hm, 'lower_upper_minus_b':hp-bmstar,
}
ineq_report=[]
for key,ex in ineqs.items():
    n,den=map(s.factor,s.fraction(s.cancel(ex)))
    assert s.factor(n.subs(q,0))==0, key
    sign=s.sign(den.subs(q,s.Rational(1,2)))
    assert sign in [1,-1], key
    assert s.Poly(den,q).count_roots(s.Rational(0),s.Rational(1))==0, key
    p=s.Poly(s.cancel(n/q)*sign,q)
    coeffs=list(reversed(p.all_coeffs()))
    lower_bound=s.factor(coeffs[0] + sum(min(s.Integer(0),x) for x in coeffs[1:]))
    assert lower_bound>0, key
    ineq_report.append({'inequality':key,'normalized_P_coefficients_ascending':[str(v) for v in coeffs],
       'P_lower_bound_0_1':str(lower_bound),'denominator_sign':str(sign),
       'denominator_roots_closed_0_1':0,'status':'PASS_FOR_ALL_0_LE_Q_LE_1'})

# Symbolic exact gap expressions match printed R5 support numbers.
assert s.factor(Dp.subs(q,1)-s.Rational(1855736,2541891)) == 0
assert s.factor(Dm.subs(q,1)-s.Rational(3695632,7625673)) == 0
for gap in (Dp,Dm):
    assert s.factor(gap.subs(q,0))==0
    assert float(gap.subs(q,s.Rational(1,2)))>0

# Independently evaluate by exact Fractions all rational cut values; no interpolation.
def FF(expr, val):
    y=s.cancel(expr.subs(q,s.Rational(val.numerator,val.denominator)))
    return Fraction(int(y.p),int(y.q))

def state(aa,bb,vs):
    det=3*aa*bb+2*aa+2*bb+1
    return [(17*aa*bb+15*aa+20*bb+10)*vs/(5*det),
            (17*aa*bb+15*aa+7*bb+5)*vs/(5*det),
            (17*aa*bb+4*aa+7*bb+2)*vs/(5*det)]

def fmt(x): return f'{x.numerator}/{x.denominator}'

rows=[]
for j in range(101):
    qc=Fraction(j,100); lo=Fraction(1)-qc/Fraction(5);hi=Fraction(1)+qc/Fraction(5)
    uv=[FF(z,qc) for z in U];lv=[FF(z,qc) for z in L]
    qplus=[Fraction(5)+qc, Fraction(5)+qc, Fraction(5)+qc]
    # complete source values are (2,1,2/5)*(5+-q), not uniform
    for case,t,src,aa,bb,gap,env in [
         ('upper',[FF(z,qc) for z in Up],[2*(5+qc),5+qc,Fraction(2,5)*(5+qc)],FF(apstar,qc),FF(bpstar,qc),FF(Dp,qc),uv),
         ('lower',[FF(z,qc) for z in Lp],[2*(5-qc),5-qc,Fraction(2,5)*(5-qc)],FF(amstar,qc),FF(bmstar,qc),FF(Dm,qc),lv)]:
        assert lo<=aa<=hi and lo<=bb<=hi, (j,case,aa,bb)
        # Exact local conservation:
        Kx=[(1+aa)*t[0]-aa*t[1], -aa*t[0]+(1+aa+bb)*t[1]-bb*t[2], -bb*t[1]+(1+bb)*t[2]]
        assert Kx==src, (j,case,'local')
        d=[abs(e-z) for e,z in zip(env,t)]
        assert all(x==gap/3 for x in d), (j,case,'equal-repair')
        assert sum(d)==gap and max(d)==gap/3
        assert sum(t)==sum(src)
        rows.append({'q':fmt(qc),'alpha':fmt(1-qc),'corner':case,
            'a_star':fmt(aa),'b_star':fmt(bb),
            'd1_K':fmt(gap),'dinf_K':fmt(gap/3),
            't1_K':fmt(t[0]),'t2_K':fmt(t[1]),'t3_K':fmt(t[2]),'status':'PASS'})

csv_path=OUT/'R6_3node_101cut_sharp_repair.csv'
with csv_path.open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
ineq_path=OUT/'R6_3node_8_inequality_certificates.csv'
with ineq_path.open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(ineq_report[0]));w.writeheader();w.writerows(ineq_report)

cert={
 'title':'Three-node exact all-alpha L1 and Linf sharp feasible-repair certificate',
 'original_source':'R5 unchanged synthetic three-node independent inter-node conductances, all ambient g_i=1 W/K',
 'q_domain':'0 <= q = 1-alpha <= 1',
 'upper_box_gap_W':str(Dp),'lower_box_gap_W':str(Dm),
 'upper_L1_K':str(Dp),'lower_L1_K':str(Dm),
 'upper_Linf_K':str(s.factor(Dp/3)),'lower_Linf_K':str(s.factor(Dm/3)),
 'upper_star_a':str(apstar),'upper_star_b':str(bpstar),
 'lower_star_a':str(amstar),'lower_star_b':str(bmstar),
 'upper_support_d1_K':str(Dp.subs(q,1)),'upper_support_dinf_K':str((Dp/3).subs(q,1)),
 'lower_support_d1_K':str(Dm.subs(q,1)),'lower_support_dinf_K':str((Dm/3).subs(q,1)),
 'upper_support_optimal_a':str(apstar.subs(q,1)),
 'upper_support_optimal_b':str(bpstar.subs(q,1)),
 'lower_support_optimal_a':str(amstar.subs(q,1)),
 'lower_support_optimal_b':str(bmstar.subs(q,1)),
 'whole_alpha_link_bounds_inequalities':len(ineq_report),
 'whole_alpha_all_8_inequalities_strict_for_q_positive':True,
 'whole_alpha_local_substitutions_verified':2,
 'independent_exact_rational_grid_checks':len(rows),
 'all_alpha_proof':'conservation hyperplane L1/Linf lower bounds; explicit a*(q),b*(q); eight polynomial positive lower certificates of link bounds; exact local substitution',
 'is_experimental_validation':False,
 'scientific_status':'exact certificate for declared synthetic R5 three-node model only; known dependency phenomenon and norm lower bounds not claimed novel universally'
}
(OUT/'R6_3node_summary.json').write_text(json.dumps(cert,indent=2,ensure_ascii=False,sort_keys=True)+'\n',encoding='utf-8')
print('PASS',len(ineq_report),'whole-domain inequalities;',len(rows),'rational reconstructed states; exact local algebra')
print('support upper d1=', cert['upper_support_d1_K'], 'dinf=',cert['upper_support_dinf_K'])
print('support lower d1=', cert['lower_support_d1_K'], 'dinf=',cert['lower_support_dinf_K'])
