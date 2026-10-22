"""Verify the additive manuscript archive without rerunning inference."""
import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
MANIFEST = Path(__file__).with_name("ADDENDUM_SHA256SUMS.txt")


def require(condition, message):
    if not condition:
        raise SystemExit("FAIL: " + message)


count = 0
for line in MANIFEST.read_text().splitlines():
    expected, name = line.split("  ", 1)
    path = ROOT / name
    require(path.resolve().is_relative_to(ROOT.resolve()), "path outside root: " + name)
    require(path.is_file(), "missing file: " + name)
    with path.open("rb") as handle:
        actual = hashlib.file_digest(handle, "sha256").hexdigest()
    require(actual == expected, "checksum mismatch: " + name)
    count += 1

fig = ROOT / "paper/figure_package/updates"
identities = fig / "Figure2_20261012/source_data/panels/Figure2D_feature_identities.csv"
with identities.open(newline="") as handle:
    rows = list(csv.DictReader(handle))
require([r["display_id"] for r in rows] == ["I15", "I06", "I19"], "Figure 2D IDs")

with (ROOT / "paper/provenance/ovarian/rctd/spacexr_DESCRIPTION").open() as handle:
    require("Version: 2.2.1" in handle.read().splitlines(), "archived spacexr version")

for base, name in [("Figure1_Atera_20261016", "Figure1"), ("Figure2_20261012", "Figure2")]:
    for suffix in ("pdf", "svg", "png"):
        path = fig / base / "figures" / (name + "." + suffix)
        require(path.is_file() and path.stat().st_size > 0, "missing figure: " + str(path))

print(f"PASS: {count} addendum checksums; main-figure files and Figure 2D identities verified")
print("This verifies archived files, not new inference, model fitting or biological validity.")
