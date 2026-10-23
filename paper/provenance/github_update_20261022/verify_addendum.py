"""The dated addendum is now integrated into the two-package distribution."""
import runpy
from pathlib import Path
root=Path(__file__).resolve().parents[3]
runpy.run_path(str(root/'tools/verify_software.py'),run_name='__main__')
