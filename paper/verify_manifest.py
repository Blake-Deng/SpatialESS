"""Current software manifest; companion objects are verified separately."""
import runpy
from pathlib import Path
runpy.run_path(str(Path(__file__).resolve().parents[1]/'tools/verify_software.py'),run_name='__main__')
