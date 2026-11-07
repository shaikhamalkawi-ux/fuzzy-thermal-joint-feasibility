#!/usr/bin/env python3
"""Exact-rational matched thermal-envelope certificate: 2 and 3 nodes.

All inputs are designed symbolic illustrations, NOT observed material data.
All-alpha monotonicity/derivative proofs are written in the report and checked
using SymPy if present; grid checks complement, never replace, those proofs.
"""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import csv, json, hashlib

HERE=Path(__file__).resolve().parent

def invert_solve(A,b):
    n=len(b)
    M=[[F(A[i][j]) for j in range(n)]+[F(b[i])] for i in range(n)]
    for j in range(n):
        k=next(i for i in range(j,n) if M[i][j]!=0)
        M[j],M[k]=M[k],M[j]
        d=M[j][j]
        M[j]=[x/d for x in M[j]]
        for i in range(n):
            if i==j: continue
            z=M[i][j]
            M[i]=[p-z*r for p,r in zip(M[i],M[j])]
    return [M[i][-1] for i in range(n)]

def evaluate(n,q, a, b=None, Q=None):
    assert n in (2,3)
    if n==2:
        assert Q is not None and len(Q)==2
        K=[[1+a,-a],[-a,1+a]]
    else:
        assert b is not None and Q is not None and len(Q)==3
        K=[[1+a,-a,F(0)],[-a,1+a+b,-b],[F(0),-b,1+b]]
    T=invert_solve(K,Q)
    assert sum(T)==sum(Q) # global energy conservation (g_i=1)
    assert all(sum(K[i][j]*T[j] for j in range(n))==Q[i] for i in range(n))
    return T

def bounds(n,q):
    lo,hi=F(1)-q/5,F(1)+q/5
    if n==2:
        ql=[F(10)-2*q,F(5)-q]; qu=[F(10)+2*q,F(5)+q]
        # T1 decreasing h, T2 increasing h, both increasing in Q
        L=[evaluate(2,q,hi,Q=ql)[0],evaluate(2,q,lo,Q=ql)[1]]
        U=[evaluate(2,q,lo,Q=qu)[0],evaluate(2,q,hi,Q=qu)[1]]
        return L,U,ql,qu
    # Q1>Q2>Q3 for all admissible Q, hence T1>T2>T3
    ql=[2*(5-q),5-q,F(2,5)*(5-q)]
    qu=[2*(5+q),5+q,F(2,5)*(5+q)]
    # exact all-alpha signs:
    # T1 decreases in a,b; T2 increases in a and decreases in b;
    # T3 increases in a,b.
    L=[evaluate(3,q,hi,hi,ql)[0],evaluate(3,q,lo,hi,ql)[1],evaluate(3,q,lo,lo,ql)[2]]
    U=[evaluate(3,q,lo,lo,qu)[0],evaluate(3,q,hi,lo,qu)[1],evaluate(3,q,hi,hi,qu)[2]]
    return L,U,ql,qu

def verify_at_cut(n,q):
    L,U,ql,qu=bounds(n,q)
    for i in range(n):
        assert L[i]<=U[i]
    # Exact, finite corner sanity checks: numerical grids do not establish all-alpha.
    ab=[F(1)-q/5,F(1)+q/5]
    hvals=list(product(ab,repeat=1 if n==2 else 2))
    qs=list(product(*[(x,y) for x,y in zip(ql,qu)]))
    corners=[evaluate(n,q,*hh,Q=list(Q)) for hh in hvals for Q in qs]
    assert len(corners)==(8 if n==2 else 32)
    for i in range(n):
        assert min(T[i] for T in corners)==L[i]
        assert max(T[i] for T in corners)==U[i]
    Qlow,Qhigh=sum(ql),sum(qu)
    gaplow=Qlow-sum(L)
    gaphigh=sum(U)-Qhigh
    if q==0:
        assert gaplow==0 and gaphigh==0
    else:
        assert gaplow>0 and gaphigh>0
    if n==2:
        exact_gaphigh=F(10)*q*(5+q)/(225-4*q*q)
        exact_gaplow=F(10)*q*(5-q)/(225-4*q*q)
        assert gaphigh==exact_gaphigh and gaplow==exact_gaplow
    return L,U,gaplow,gaphigh,len(corners)

def symbolic_proof():
    from sympy import symbols, Matrix, Rational, diff, factor, simplify
    q,a,b=symbols('q a b',real=True)
    K=Matrix([[1+a,-a,0],[-a,1+a+b,-b],[0,-b,1+b]])
    v=Matrix([2,1,Rational(2,5)])
    T=K.inv()*v
    D=3*a*b+2*a+2*b+1
    assert factor(K.det()-D)==0
    # sign factorisations: each multiplier is strictly positive at a,b>0,
    # so signs of every partial derivative follow from the exact identities.
    factors=[
      [-(2*b+1)*(13*b+5),-a*(11*a+3)],
      [(b+1)*(13*b+5),-(a+1)*(11*a+3)],
      [b*(13*b+5),(2*a+1)*(11*a+3)]
    ]
    for i in range(3):
      for j,p in enumerate((a,b)):
        assert simplify(diff(T[i],p)-factors[i][j]/(5*D*D))==0
    # vector Q is (5 +/- q)*v: signs persist over all levels q in [0,1]
    Hlo=1-q/5;Hhi=1+q/5
    lower=[(Hhi,Hhi),(Hlo,Hhi),(Hlo,Hlo)]
    upper=[(Hlo,Hlo),(Hhi,Hlo),(Hhi,Hhi)]
    TsL=[factor((5-q)*T[i].subs({a:x,b:y},simultaneous=True)) for i,(x,y) in enumerate(lower)]
    TsU=[factor((5+q)*T[i].subs({a:x,b:y},simultaneous=True)) for i,(x,y) in enumerate(upper)]
    for i in range(3):
        assert factor(TsL[i].subs(q,0)-TsU[i].subs(q,0))==0
        assert TsU[i].subs(q,1)>TsL[i].subs(q,1)
    dlow=factor(Rational(17,5)*(5-q)-sum(TsL))
    dhi=factor(sum(TsU)-Rational(17,5)*(5+q))
    assert dlow.subs(q,1)==Rational(3695632,7625673)
    assert dhi.subs(q,1)==Rational(1855736,2541891)
    return {
      '3node_det':str(D),
      '3node_symbolic_derivative_signs':[['negative','negative'],['positive','negative'],['positive','positive']],
      '3node_global_lower_gap_q':str(dlow),
      '3node_global_upper_gap_q':str(dhi),
      '3node_bounds':{f'T{i+1}_lower(q)':str(TsL[i]) for i in range(3)}|{f'T{i+1}_upper(q)':str(TsU[i]) for i in range(3)},
      'symbolic_tests_passed': '1 determinant + 6 signed derivatives + 6 boundary/width tests + 2 support gaps',
    }

def main():
    rows=[]; checks=0
    for n in (2,3):
        prev=None
        for k in range(101):
            alpha=F(k,100);q=1-alpha
            L,U,dl,du,cc=verify_at_cut(n,q)
            checks+=cc
            if prev:
                pL,pU=prev
                assert all(pL[i]<=L[i]<=U[i]<=pU[i] for i in range(n))
            prev=L,U
            rows.append({
                'nodes':str(n),'alpha':str(alpha),'q':str(q),
                **{f'T{i+1}_min_K':str(L[i]) for i in range(n)},
                **{f'T{i+1}_max_K':str(U[i]) for i in range(n)},
                'source_lower_W':str(sum(bounds(n,q)[2])),
                'source_upper_W':str(sum(bounds(n,q)[3])),
                'energy_deficit_low_corner_W':str(dl),
                'energy_excess_high_corner_W':str(du),
                'independent_corner_cases_checked':str(cc),
            })
    cols=['nodes','alpha','q','T1_min_K','T1_max_K','T2_min_K','T2_max_K','T3_min_K','T3_max_K',
      'source_lower_W','source_upper_W','energy_deficit_low_corner_W','energy_excess_high_corner_W','independent_corner_cases_checked']
    out=HERE/'thermal_two_three_node_101cuts.csv'
    with out.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=cols);w.writeheader();w.writerows(rows)
    theory=symbolic_proof()
    data={
       'research_status':'illustrative mathematical certificate, NOT measured or validated material performance',
       'certified_claim':'The componentwise envelope contains energy-infeasible corners for every alpha<1 (both networks).',
       '2node_upper_gap': '10*q*(5+q)/(225-4*q*q)',
       '2node_lower_gap': '10*q*(5-q)/(225-4*q*q)',
       '2node_support_deficit_W':'40/221','2node_support_excess_W':'60/221',
       '3node_support_deficit_W':'3695632/7625673',
       '3node_support_excess_W':'1855736/2541891',
       '101_cut_checks':202,
       'exact_conserved_corner_solutions_checked':checks,
       'monotone_nested_cut_transitions_checked':200,
       **theory,
    }
    (HERE/'symbolic_certificate.json').write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n')
    print('PASS',json.dumps({k:data[k] for k in ['101_cut_checks','exact_conserved_corner_solutions_checked','monotone_nested_cut_transitions_checked']}))
    print('Two-node q=1 gap:',data['2node_support_deficit_W'],data['2node_support_excess_W'])
    print('Three-node q=1 gap:',data['3node_support_deficit_W'],data['3node_support_excess_W'])

if __name__=='__main__': main()
