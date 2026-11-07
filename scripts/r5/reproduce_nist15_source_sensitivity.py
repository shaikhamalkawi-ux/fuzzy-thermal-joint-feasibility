#!/usr/bin/env python3
"""NIST official published historical table: source consistency and descriptive sensitivity ONLY.
Original NIST SRM 1450d certificate: DOI 10.6028/NIST.SP.260-173.
Published 15 original results: NIST JRES (2014) Table I1, PDF page 65 (zero-index 64),
DOI 10.6028/jres.119.012; official source datasets ZIP advertised at
https://www.nist.gov/el/beed/nist-srm-1450 but binary not acquired in this package.
Not holdout, not a recalibration, not validation of a thermal network, and not a statistical
confidence or fuzzy membership inference.
"""
from fractions import Fraction as F
from pathlib import Path
from collections import defaultdict
import csv, json, math
BASE=Path(__file__).resolve().parent
A=F('0.000110489')  # NIST certified mean-temperature conductivity fit
G0=F(331467,1016000) # fixed hypothetical A/L from R4
SIG=F(1,5)
ETA_ILLUSTRATIVE=F(1,100)

def dec(x): return float(x)
def gaps(e):
    plus=G0*(1+SIG)*e/(3*(3+e))
    minus=G0*(1-SIG)*e/(3*(3-e))
    dp=(1+SIG)*e/(3*(1-e)*(3+e))
    dm=(1-SIG)*e/(3*(1+e)*(3-e))
    return {"upper_box_excess_W":plus,"lower_box_deficit_W":minus,
            "upper_sharp_L1_potential_K":dp,"upper_sharp_Linfty_potential_K":dp/2,
            "lower_sharp_L1_potential_K":dm,"lower_sharp_Linfty_potential_K":dm/2}

rows=[]; bygrp=defaultdict(list)
with (BASE/'NIST_1450d_TableI1_source_rows.csv').open(newline='',encoding='utf8') as f:
    for r in csv.DictReader(f):
        lam,Tm=F(r['lambda_exp_W_mK']),F(r['Tmean_K']); pred=A*Tm
        residual=lam/pred-1
        row=dict(r)
        row.update(lambda_pred_W_mK=f'{float(pred):.10f}',
            relative_residual_pct=f'{float(100*residual):.9f}',
            absolute_relative_residual_pct=f'{float(abs(100*residual)):.9f}',
            signed_residual_exact=f'{residual.numerator}/{residual.denominator}')
        rows.append((row,residual));bygrp[r['nominal_T_group_K']].append((r,residual))
assert len(rows)==15 and sorted(map(int,bygrp))==[280,295,310,325,340]
assert all(len(v)==3 for v in bygrp.values())
argmax, eta_signed=max(rows,key=lambda rr:abs(rr[1]))
eta=abs(eta_signed); assert eta < ETA_ILLUSTRATIVE
res=[d for row,d in rows]
mean=sum(res)/15
mae=sum(map(abs,res))/15
rmse=math.sqrt(sum(float(x*x) for x in res)/15)
within=sum(abs(x)<=F(1,100) for x in res)
assert within == 15
assert sum(F(r['Tmean_K'])<280 for r,d in rows) == 2
# source check: first and last row transcribed exactly from publisher Table I1
assert rows[0][0]['lambda_exp_W_mK']=='0.031123'
assert rows[-1][0]['lambda_exp_W_mK']=='0.037817'

variants=[]
for label,e in [('nominal_zero_variation',F(0)),('max_observed_abs_residual_descriptive',eta),('R4_illustrative_one_percent',ETA_ILLUSTRATIVE)]:
    out={'scenario':label,'eta_exact':f'{e.numerator}/{e.denominator}', 'eta_decimal':float(e)}
    for key,v in gaps(e).items():
        out[key]=f'{v.numerator}/{v.denominator}'
        out[key+'_decimal']=float(v)
    variants.append(out)
with (BASE/'NIST_1450d_TableI1_residuals_recomputed.csv').open('w',newline='',encoding='utf8') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0][0]));w.writeheader();w.writerows([r for r,d in rows])
with (BASE/'NIST_sensitivity_results.csv').open('w',newline='',encoding='utf8') as f:
    w=csv.DictWriter(f,fieldnames=list(variants[0]));w.writeheader();w.writerows(variants)

grouped=[]
for g in sorted(bygrp,key=int):
    data=bygrp[g]
    grouped.append({'nominal_temperature_K':int(g),'n':len(data),
        'min_observed_lambda_W_mK':min(float(r['lambda_exp_W_mK']) for r,d in data),
        'max_observed_lambda_W_mK':max(float(r['lambda_exp_W_mK']) for r,d in data),
        'max_abs_residual_pct':max(float(abs(100*d)) for r,d in data)})

cert={'source':'Zarr, Heckert, Leigh 2014, DOI:10.6028/jres.119.012, Table I1, SRM1450d, Lot2009',
    'certification_source':'NIST SP260-173 DOI:10.6028/NIST.SP.260-173',
    'input_access':'published source Table I1 transcribed from official NIST PDF; NIST source-data ZIP binary NOT obtained',
    'data_points':15,'number_temperature_groups':5,
    'count_nominally_below_280_K_due_to_printed_precision':2,
    'all_abs_deviations_within_1_percent':within,
    'signed_mean_percent':float(100*mean),
    'mean_abs_deviation_percent':float(100*mae),
    'rms_relative_deviation_percent':100*rmse,
    'max_abs_deviation_percent':float(100*eta),
    'max_abs_deviation_row':int(argmax['row']),
    'max_abs_deviation_exact_fraction':f'{eta.numerator}/{eta.denominator}',
    'source_temperature_group_comparison':grouped,
    'matched_marginal_sensitivities':variants,
    'guardrails':["These same-lot historical observations partly underlie SRM1450d certification; NOT independent validation.",
                 "The 1% NIST expanded uncertainty with coverage factor ~2 is NOT a hard deterministic bound or a fuzzy membership function.",
                 "The maximum observed absolute residual is a finite-sample descriptive statistic, not an inferential bound.",
                 "The 1% and empirical observed residual bands are hypothetical model sensitivity inputs, not calibrated uncertainty distributions.",
                 "Path areas, lengths, process loads, multi-node topology and local Fourier law remain hypothetical/conditional.",
                 "Temperature validity 280-340 K is for mean-temperature relation; two raw printed mean temperatures are 279.998 K.",
                 "Every real network state in R4 conserves energy; positive gaps belong to unattainable recombinations of marginal coordinates."]}
(BASE/'NIST15_source_audit_certificate.json').write_text(json.dumps(cert,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
print(json.dumps({k:v for k,v in cert.items() if k in ['data_points','signed_mean_percent','mean_abs_deviation_percent','rms_relative_deviation_percent','max_abs_deviation_percent','max_abs_deviation_row','all_abs_deviations_within_1_percent']},indent=2))
for v in variants:
 print(v['scenario'],'eta',v['eta_decimal'],'gap+',v['upper_box_excess_W_decimal'],'gap-',v['lower_box_deficit_W_decimal'])
