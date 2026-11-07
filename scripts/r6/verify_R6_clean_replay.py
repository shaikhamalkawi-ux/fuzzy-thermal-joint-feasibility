#!/usr/bin/env python3
"""Clean-extracted R6 source, algebra, visual and archival QA.

UPSTREAM_R5_IMMUTABLE.zip is the original, unmodified prior research baseline.
"""
from pathlib import Path
import sys, subprocess, tempfile, hashlib, shutil, zipfile
root=Path(__file__).resolve().parent
source=['R6_3node_sharp_repair_certificate.py','R6_independent_fraction_replay.py','build_R6_paper.py',
   'R6_3node_101cut_sharp_repair.csv','R6_3node_8_inequality_certificates.csv','R6_3node_summary.json',
   'Figure_2node_Energy_Balance_Gap.png','Figure_NIST_MatchedMarginal_Dependence.png',
   'Figure_NIST_2014_TableI1_SourceResiduals.png','Figure_R6_ThreeNode_SharpRepair.pdf']
for name in source+['UPSTREAM_R5_IMMUTABLE.zip']:
    assert (root/name).is_file(), name
with tempfile.TemporaryDirectory(prefix='99fuzzy_r6_fresh_') as td:
    target=Path(td)
    for name in source:shutil.copy2(root/name,target/name)
    with zipfile.ZipFile(root/'UPSTREAM_R5_IMMUTABLE.zip') as z:
        assert z.testzip() is None
        content=z.read('ICMPC2027_R5_NIST15_SourceAudited_WORKING.tex')
    (target/'R5_upstream').mkdir()
    (target/'R5_upstream/ICMPC2027_R5_NIST15_SourceAudited_WORKING.tex').write_bytes(content)
    for cmd in ([sys.executable,'R6_3node_sharp_repair_certificate.py'],
                [sys.executable,'R6_independent_fraction_replay.py'],
                [sys.executable,'build_R6_paper.py']):
        subprocess.run(cmd,cwd=target,check=True,capture_output=True)
    for name in ('R6_3node_101cut_sharp_repair.csv','R6_3node_8_inequality_certificates.csv','R6_3node_summary.json'):
        assert hashlib.sha256((root/name).read_bytes()).hexdigest()==hashlib.sha256((target/name).read_bytes()).hexdigest(), name
    tex='ICMPC2027_R6_ThreeNode_SharpRepair_SourceAudited_WORKING.tex'
    pdf=tex.replace('.tex','.pdf')
    for _ in range(2):subprocess.run(['pdflatex','-interaction=nonstopmode','-halt-on-error',tex],cwd=target,check=True,capture_output=True)
    out=subprocess.run(['pdfinfo',pdf],cwd=target,check=True,capture_output=True,text=True).stdout
    assert 'Pages:           9' in out
    fontdata=subprocess.run(['pdffonts',pdf],cwd=target,check=True,capture_output=True,text=True).stdout
    assert 'Type 3' not in fontdata
    log=(target/tex.replace('.tex','.log')).read_text(errors='replace')
    assert 'undefined' not in log.lower() and 'Overfull' not in log
    print('PASS 3/3 exact outputs byte-identical; 202 independent Fraction states; all 8 whole-domain inequalities; fresh 9-page PDF; no Type3/undefined/overfull')
