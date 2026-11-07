#!/usr/bin/env python3
"""Public source-level replay of 19 original R4/R5/R6/G2 generated outputs.

Always runs scripts from clean temporary working directories. The exact
source-defined checks are mathematical/archival, not field validation.
"""
from __future__ import annotations
from pathlib import Path
from tempfile import TemporaryDirectory
from shutil import copy2
import subprocess
import hashlib
import json
import sys

ROOT=Path(__file__).resolve().parent

OUTPUTS={
 'r4':['R4_exact_sharp_repair_101cuts.csv','R4_exact_sharp_repair_certificate.json',
       'thermal_two_three_node_101cuts.csv','symbolic_certificate.json',
       'NIST_commonK_101cut_certificate.csv','NIST_commonK_certificate.json',
       'matched_marginal_101cuts.csv','matched_marginal_certificate.json',
       'matched_marginal_sympy_certificate.json'],
 'r5':['NIST_1450d_TableI1_residuals_recomputed.csv','NIST_sensitivity_results.csv',
       'NIST15_source_audit_certificate.json','NIST15_descriptive_eta_101cut_grid.csv',
       'NIST15_independent_corner_certificate.json'],
 'r6':['R6_3node_101cut_sharp_repair.csv','R6_3node_8_inequality_certificates.csv',
       'R6_3node_summary.json'],
 'g2':['outputs/G2_INDEPENDENT_EIGHT_CERTIFICATES.csv',
       'outputs/G2_INDEPENDENT_EIGHT_CERTIFICATES.json']
}
SCRIPTS={
 'r4':['r4_exact_sharp_repair.py','reproduce_exact_thermal.py','reproduce_nist_common_k.py',
       'reproduce_matched_marginal_dependence.py','independent_symbolic_proof.py'],
 'r5':['reproduce_nist15_source_sensitivity.py','independent_3232_cut_checker.py'],
 'r6':['R6_3node_sharp_repair_certificate.py','R6_independent_fraction_replay.py'],
 'g2':['scripts/independent_certificate_audit.py']
}

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def validate_source_hashes():
    record=json.loads((ROOT/'provenance/ORIGINAL_FILE_SHA256.json').read_text())['files']
    for name,metadata in record.items():
        p=ROOT/name
        assert p.is_file() and digest(p)==metadata['sha256'], 'source mismatch: '+name
    return len(record)

def reproduce():
    checked=[]; log=[]
    for release in OUTPUTS:
        with TemporaryDirectory(prefix='99fuzzyc_'+release+'_') as tmp:
            work=Path(tmp)
            src=ROOT/'scripts'/release
            expected=ROOT/'data'/release
            # Copy only source, not generated outputs. The NIST historical
            # transcribed table is the sole external input of the R5 scripts.
            if release!='g2':
                for script in src.glob('*.py'):copy2(script,work/script.name)
                if release=='r5':
                    copy2(expected/'NIST_1450d_TableI1_source_rows.csv',
                          work/'NIST_1450d_TableI1_source_rows.csv')
            else:
                (work/'scripts').mkdir()
                copy2(src/'independent_certificate_audit.py',work/'scripts/independent_certificate_audit.py')
                (work/'inputs/r6').mkdir(parents=True)
                for input_file in (expected/'inputs/r6').iterdir():
                    copy2(input_file,work/'inputs/r6'/input_file.name)
                (work/'outputs').mkdir()
            for script in SCRIPTS[release]:
                command=[sys.executable,script]
                result=subprocess.run(command,cwd=work,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,timeout=180,env={'PYTHONDONTWRITEBYTECODE':'1','PYTHONPATH':''})
                log.append({'release':release,'script':script,'status':result.returncode,
                            'stdout_tail':result.stdout[-700:],'stderr_tail':result.stderr[-300:]})
                if result.returncode:
                    raise RuntimeError(f'{release}: {script} failed: {result.stderr[-1800:]} {result.stdout[-1800:]}')
            for output in OUTPUTS[release]:
                actual=work/output;baseline=expected/output
                if not actual.is_file():raise AssertionError('missing regenerated output '+str(actual))
                if not baseline.is_file():raise AssertionError('missing source benchmark '+str(baseline))
                if digest(actual)!=digest(baseline):
                    raise AssertionError('non-identical output: '+release+'/'+output)
                checked.append(release+'/'+output)
    assert len(checked)==19
    return checked,log

if __name__=='__main__':
    try:
        source_count=validate_source_hashes()
        results,log=reproduce()
        print(json.dumps({'result':'PASS','source_hash_checks':source_count,
            'byte_identical_outputs':len(results),
            'groups':{key:len(OUTPUTS[key]) for key in OUTPUTS},
            'claim_boundary':'Synthetic, exact stated thermal-network models; no experimental validation',
            'outputs':results},indent=2))
    except Exception as exc:
        print('FAIL:',exc,file=sys.stderr)
        sys.exit(1)
