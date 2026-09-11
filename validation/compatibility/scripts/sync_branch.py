"""Copy the tested candidate engine into the standalone bridge repository."""
import shutil
from pathlib import Path

root = Path(__file__).resolve().parent
src = root / 'work/SpatialESS'
dst = root / 'work/SpatialESS-SPARKLE/vendor/SpatialESS'
assert dst.resolve().is_relative_to((root / 'work').resolve())
shutil.rmtree(dst)
shutil.copytree(src, dst, ignore=shutil.ignore_patterns('*.o', '*.so', '*.dll', '.git', '__pycache__', '*.pyc'))
helper = Path('/home/dzf/cellchat_acceleration/validation/collaborator_cervical_20260911/export_collaborator_debug.R')
for repo in (src, dst):
    shutil.copy2(helper, repo / 'benchmarks/export_spatialcellchat_diagnostic.R')
print('Bundled core synchronized; no external main checkout required.')
