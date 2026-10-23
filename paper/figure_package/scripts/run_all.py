"""Current saved-result entry points. Mount the companion first or set SPATIALESS_MATERIALS."""
from pathlib import Path
import subprocess,sys,os
root=Path(__file__).resolve().parents[3]
for name in ['Figure1','Figure2','S1','S2','S3','S4']:
    cmd=[sys.executable,str(root/'tools/redraw.py'),name]
    if os.environ.get('SPATIALESS_MATERIALS'):cmd += ['--materials',os.environ['SPATIALESS_MATERIALS']]
    subprocess.run(cmd,cwd=root,check=True)
