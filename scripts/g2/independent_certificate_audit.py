#!/usr/bin/env python3
"""99FUZZYC — independent targeted symbolic check, not human validation.

Derives the 3-node link repair entirely from the original balance matrix and
source vector, not by importing the R6 solver or its certificate table. It
compares to original archived certificates *after* independently deriving them.
All numbers are exact SymPy rationals. No grid is used to prove continuum signs.
"""
import csv, json, hashlib, platform
from pathlib import Path
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
ORIG=ROOT/'inputs/r6'
OUT=ROOT/'outputs'
OUT.mkdir(parents=True,exist_ok=True)

def sha256(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()

q=sp.symbols('q', real=True)
a,b=sp.symbols('a b',positive=True)
H=sp.Matrix([[1+a,-a,0],[-a,1+a+b,-b],[0,-b,1+b]])
v=sp.Matrix([2,1,sp.Rational(2,5)])
# Derive F from H instead of transcribing the published rational expressions.
F=[sp.factor(z) for z in H.LUsolve(v)]
derivatives={f'F{i+1}_d{x}':str(sp.factor(sp.diff(f,x))) for i,f in enumerate(F) for x in (a,b)}
expected_signs=[(-1,-1),(1,-1),(1,1)]
# Certificate of derivative signs: products of positive affine terms for all a,b>0.
for i,(f,signs) in enumerate(zip(F,expected_signs)):
 for x,sign in zip((a,b),signs):
  ex=sp.factor(sp.diff(f,x))
  n,d=sp.fraction(ex)
  # Multiplication by sign yields positive products for a,b>0; test polynomial factors.
  num=sp.factor(sign*n)
  assert num.is_Mul or num.is_Add or num.is_Pow or num.is_Symbol or num.is_Number
  assert sp.Poly(num,a,b).coeffs() and all(coef>0 for coef in sp.Poly(num,a,b).coeffs()), (i,x,num)
  assert all(coef>0 for coef in sp.Poly(sp.expand(d),a,b).coeffs()), (i,x,d)

lo,hi=1-q/5,1+q/5
# Argmin/argmax conductance corners follow from independently certified partial signs.
lo_pairs=[(hi,hi),(lo,hi),(lo,lo)]
hi_pairs=[(lo,lo),(hi,lo),(hi,hi)]
L=sp.Matrix([sp.factor((5-q)*F[i].subs({a:xy[0],b:xy[1]})) for i,xy in enumerate(lo_pairs)])
U=sp.Matrix([sp.factor((5+q)*F[i].subs({a:xy[0],b:xy[1]})) for i,xy in enumerate(hi_pairs)])
Qlo=(5-q)*v;Qhi=(5+q)*v
Glo=sp.factor(sum(Qlo)-sum(L))
Ghi=sp.factor(sum(U)-sum(Qhi))
# Exact strict gap positivity from factored rational expressions over the full domain.
# The common denominator (q-10)(q+10)(3q-20)(3q+20)(3q^2-200)
# has three strictly negative factors on [0,1]; hence sign -1.
Dgap=(q-10)*(q+10)*(3*q-20)*(3*q+20)*(3*q*q-200)
assert sp.Poly(Dgap,q).count_roots(0,1)==0
assert sp.sign(Dgap.subs(q,sp.Rational(1,2)))==-1
upper_P=72*q**4-5*q**3-8600*q**2+500*q+240000
lower_P=72*q**4+5*q**3-8600*q**2-500*q+240000
def coefficient_positive_lower(expr):
 cs=list(reversed(sp.Poly(expr,q).all_coeffs()))
 return cs[0]+sum(min(sp.Integer(0),x) for x in cs[1:])
assert coefficient_positive_lower(upper_P)==sp.Integer(231395)
assert coefficient_positive_lower(lower_P)==sp.Integer(230900)
assert sp.factor(Ghi+4*q*(q+5)*upper_P/Dgap)==0
assert sp.factor(Glo-4*q*(q-5)*lower_P/Dgap)==0
# Thus Ghi>0 and Glo>0 for 0<q<=1; both vanish exactly at q=0.
assert sp.factor(Ghi.subs(q,0))==0 and sp.factor(Glo.subs(q,0))==0

# Balanced candidate: equality in necessary Linf bound.
Tlo=sp.Matrix([sp.factor(z+Glo/3) for z in L]); Thi=sp.Matrix([sp.factor(z-Ghi/3) for z in U])

def reconstruct(T,Q):
 x=sp.factor((Q[0]-T[0])/(T[0]-T[1]))
 y=sp.factor((T[2]-Q[2])/(T[1]-T[2]))
 return x,y

ahi,bhi=reconstruct(Thi,Qhi)
alo,blo=reconstruct(Tlo,Qlo)

identities={}
for label,T,Q,x,y,g,env in [('upper',Thi,Qhi,ahi,bhi,Ghi,U),('lower',Tlo,Qlo,alo,blo,Glo,L)]:
 K=H.subs({a:x,b:y})
 res=[sp.factor(z) for z in K*T-Q]
 assert res==[0,0,0], (label,res)
 assert sp.factor(sum(T)-sum(Q))==0
 abs_diffs=[sp.factor(e-t) if label=='upper' else sp.factor(t-e) for e,t in zip(env,T)]
 assert [sp.factor(e-g/3) for e in abs_diffs]==[0,0,0]
 identities[label]={'K_times_T_minus_Q_exact':[str(z) for z in res],
   'sum_T_minus_sum_Q': '0',
   'three_equal_coordinate_corrections_K':[str(z) for z in abs_diffs],
   'gap_W_numeric_for_gref_1W_per_K':str(g),
   'a_star_over_gref':str(x),'b_star_over_gref':str(y)}

# The eight signed difference expressions independently derived above.
differences={
 'upper_a_minus_lower':ahi-lo,'upper_upper_minus_a':hi-ahi,
 'upper_b_minus_lower':bhi-lo,'upper_upper_minus_b':hi-bhi,
 'lower_a_minus_lower':alo-lo,'lower_upper_minus_a':hi-alo,
 'lower_b_minus_lower':blo-lo,'lower_upper_minus_b':hi-blo}

# Lower-bound positivity certificate:
# when P(q)=sum_k c_k q^k, 0<=q<=1, P>=c_0+sum_{k>=1}min(0,c_k).
# Negative terms are bounded below by c_k; positive terms bounded below by 0.
recs=[]
for name,ex in differences.items():
 n,d=sp.fraction(sp.cancel(ex));n=sp.factor(n);d=sp.factor(d)
 assert sp.factor(n.subs(q,0))==0, name
 dpoly=sp.Poly(d,q)
 assert dpoly.degree()>=0
 roots=dpoly.count_roots(0,1)
 assert roots==0, (name,roots)
 dsign=int(sp.sign(d.subs(q,sp.Rational(1,2))))
 assert dsign in (-1,1)
 normalized=sp.Poly(sp.cancel(n/q)*dsign,q)
 c=list(reversed(normalized.all_coeffs()))
 bound=sp.simplify(c[0]+sum(min(sp.Integer(0),j) for j in c[1:]))
 assert bound>0, (name,bound)
 # Need positive numerator: sign(d)*n/q=P; therefore ex = q*P/(sign(d)*d).
 assert sp.cancel(ex-q*normalized.as_expr()/(dsign*d))==0, name
 recs.append({'id':name,'original_difference':str(ex),'factored_numerator':str(n),
  'factored_denominator':str(d),'denominator_sign':dsign,
  'denominator_roots_in_closed_0_1':int(roots),
  'normalized_positive_P_coefficients_ascending':[str(j) for j in c],
  'strict_positive_lower_bound':[str(bound), 'P(q) >= P(0) + sum_{k>=1} min(0,c_k)'],
  'domain':'0 <= q <= 1', 'zero_at_q0':True, 'strict_positive_for_q_gt_0':True})

# Check all 8 against archived R6 values, without relying on them to compute.
with (ORIG/'R6_3node_8_inequality_certificates.csv').open(newline='') as f:
 original={r['inequality']:r for r in csv.DictReader(f)}
for x in recs:
 y=original[x['id']]
 assert x['normalized_positive_P_coefficients_ascending']==sp.sympify(y['normalized_P_coefficients_ascending']), x['id']
 assert x['strict_positive_lower_bound'][0]==y['P_lower_bound_0_1'], x['id']
 assert x['denominator_sign']==int(y['denominator_sign']), x['id']

with (ORIG/'R6_3node_summary.json').open() as f:summary=json.load(f)
assert str(Ghi.subs(q,1))==summary['upper_support_d1_K']
assert str(Glo.subs(q,1))==summary['lower_support_d1_K']
assert str((Ghi/3).subs(q,1))==summary['upper_support_dinf_K']
assert str((Glo/3).subs(q,1))==summary['lower_support_dinf_K']
for e,name in [(ahi,'upper_support_optimal_a'),(bhi,'upper_support_optimal_b'),(alo,'lower_support_optimal_a'),(blo,'lower_support_optimal_b')]:
 assert str(sp.cancel(e.subs(q,1)))==summary[name]

proof={
 'status':'PASS — target reconstructed from original thermal balance, NOT an external human audit',
 'project':'99FUZZYC','scientific_baseline':'R6', 'editorial_baseline':'R6P5',
 'proof_domain':'0 <= q=1-alpha <=1, 4/5 <= a,b <=6/5; node conductances g_ref=1 W/K',
 'matrix_H':[[str(z) for z in row] for row in H.tolist()],
 'source_vector_v':[str(z) for z in v],
 'F_from_symbolic_LU_solve':[str(z) for z in F],
 'six_sign_factor_derivatives':derivatives,
 'whole_domain_strict_gap_positivity':{'common_denominator':str(sp.factor(Dgap)),'denominator_sign':-1, 'upper_P':str(upper_P), 'upper_P_lower_bound':str(coefficient_positive_lower(upper_P)), 'lower_P':str(lower_P), 'lower_P_lower_bound':str(coefficient_positive_lower(lower_P)), 'conclusion':'G_upper(q)>0 and G_lower(q)>0 for 0<q<=1 and both equal zero at q=0'},
 'upper_marginal_coordinates':[str(z) for z in U],
 'lower_marginal_coordinates':[str(z) for z in L],
 'identities':identities,
 'eight_whole_domain_link_inequalities':recs,
 'support_gaps_W':{'upper':str(Ghi.subs(q,1)),'lower':str(Glo.subs(q,1))},
 'proof_method':'Exact symbolic identities and integer coefficient lower bounds; denominator root count exact; finite sampling not used to justify all q',
 'archive_reference':{
  'original_r6_zip_sha256':'2f1583a67f106c6ea70ab72140c56133afa7a1569f3572a09dae1c40d84a192a',
  'r6_csv_sha256':sha256(ORIG/'R6_3node_8_inequality_certificates.csv'),
  'r6_script_sha256':sha256(ORIG/'R6_3node_sharp_repair_certificate.py')},
 'limitations':'Case-specific 2/3-node synthetic thermal model only; no empirical validation, no external human proof adjudication, no universal statement',
 'environment':{'python':platform.python_version(),'sympy':sp.__version__}
}
(OUT/'G2_INDEPENDENT_EIGHT_CERTIFICATES.json').write_text(json.dumps(proof,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
with (OUT/'G2_INDEPENDENT_EIGHT_CERTIFICATES.csv').open('w',newline='',encoding='utf-8') as f:
 w=csv.DictWriter(f,['id','denominator_sign','denominator_roots_in_closed_0_1','normalized_positive_P_coefficients_ascending','strict_positive_lower_bound','factored_numerator','factored_denominator']);w.writeheader();
 for x in recs:w.writerow({k:json.dumps(x[k]) if isinstance(x[k],list) else x[k] for k in w.fieldnames})
print('G2 INDEPENDENT TARGETED CERTIFICATE PASS:',len(recs),'inequalities; 6 derivative signs; 2 exact energy identities; 2 local balance identities')
for r in recs: print(r['id'], 'P_lower_bound=', r['strict_positive_lower_bound'][0], 'den_sign=',r['denominator_sign'])
print('Support gaps [W]:',proof['support_gaps_W'])
