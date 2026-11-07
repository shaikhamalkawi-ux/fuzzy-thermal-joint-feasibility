#!/usr/bin/env python3
"""Second implementation, stdlib-only independent audit of R6 CSV against source heat-balance equations.
Recreates extrema from the original interval corner definitions, without SymPy expressions.
"""
from pathlib import Path
import csv, json
from fractions import Fraction as F
root=Path(__file__).resolve().parent
path=root/'R6_3node_101cut_sharp_repair.csv'

def fraction(x):
    return F(x)

def original_projection(a,b,amp):
    d=3*a*b+2*a+2*b+1
    return [amp*(17*a*b+15*a+20*b+10)/(5*d),
            amp*(17*a*b+15*a+7*b+5)/(5*d),
            amp*(17*a*b+4*a+7*b+2)/(5*d)]

def matrix_apply(a,b,t):
    x,y,z=t
    return ((1+a)*x-a*y,-a*x+(1+a+b)*y-b*z,-b*y+(1+b)*z)
checks=0;counts={'upper':0,'lower':0}; seen=set()
with path.open(encoding='utf-8',newline='') as f:
 for row in csv.DictReader(f):
    q=fraction(row['q']);alpha=fraction(row['alpha']); assert q+alpha==1
    k=row['corner'];counts[k]+=1
    assert (q,k) not in seen;seen.add((q,k))
    a=fraction(row['a_star']);b=fraction(row['b_star']);t=tuple(fraction(row[f't{i}_K']) for i in (1,2,3))
    lm=1-q/5;um=1+q/5
    assert lm<=a<=um and lm<=b<=um
    if k=='upper':
        env=(original_projection(lm,lm,5+q)[0],original_projection(um,lm,5+q)[1],original_projection(um,um,5+q)[2])
        src=(2*(5+q),5+q,F(2,5)*(5+q))
    else:
        env=(original_projection(um,um,5-q)[0],original_projection(lm,um,5-q)[1],original_projection(lm,lm,5-q)[2])
        src=(2*(5-q),5-q,F(2,5)*(5-q))
    assert matrix_apply(a,b,t)==src
    g=sum(abs(env[i]-t[i]) for i in range(3))
    assert g==fraction(row['d1_K'])
    assert max(abs(env[i]-t[i]) for i in range(3))==fraction(row['dinf_K'])==g/3
    assert sum(t)==sum(src)
    checks+=1
assert counts=={'upper':101,'lower':101} and checks==202
print('PASS',checks,'fully independent Fraction-only exact source/model checks')
