"""Clean local replay verifies determinism of the 5 R5 numerical output files."""
from pathlib import Path
from tempfile import TemporaryDirectory
from shutil import copyfile
import subprocess, hashlib, json
SRC=Path(__file__).resolve().parent
inputs=['NIST_1450d_TableI1_source_rows.csv','reproduce_nist15_source_sensitivity.py','independent_3232_cut_checker.py']
outputs=['NIST_1450d_TableI1_residuals_recomputed.csv','NIST_sensitivity_results.csv','NIST15_source_audit_certificate.json','NIST15_descriptive_eta_101cut_grid.csv','NIST15_independent_corner_certificate.json']
with TemporaryDirectory(prefix='R5_clean_replay_') as tmp:
 p=Path(tmp)
 for f in inputs:copyfile(SRC/f,p/f)
 for script in ['reproduce_nist15_source_sensitivity.py','independent_3232_cut_checker.py']:
  subprocess.run(['python3',str(p/script)],cwd=p,check=True,capture_output=True,text=True)
 assert all((p/f).read_bytes()==(SRC/f).read_bytes() for f in outputs)
 print('R5 clean replay PASS:',len(outputs),'exact outputs byte identical')
 report={'status':'PASS','replayed_output_files':len(outputs),'files':{f:hashlib.sha256((SRC/f).read_bytes()).hexdigest() for f in outputs}}
 (SRC/'R5_CleanReplay_Checksums.json').write_text(json.dumps(report,indent=2)+'\n')
