from pathlib import Path
import re

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
full_attempt = pd.read_csv(ROOT/'source_data/current/cervical_700k_reference_attempt.tsv', sep='\t').iloc[0]
status = pd.read_csv(ROOT/'source_data/panels/Figure1I_completion.csv')
extension = status.query("tissue == 'cervical' and N in [200000, 400000]")
if len(extension) != 4:
    raise SystemExit('Expected four cervical 200k/400k status rows')

def state(n, method):
    rows = extension[(extension.N == n) & (extension.method == method)]
    if len(rows) != 1:
        raise SystemExit(f'Missing unique status for cervical {n} {method}')
    return str(rows.iloc[0].status)

ess_200 = state(200000, 'SpatialESS')
ess_400 = state(400000, 'SpatialESS')
ref_200 = state(200000, 'SpatialCellChat')
ref_400 = state(400000, 'SpatialCellChat')
state_labels = {'allocation_failure_under_memory_budget': 'allocation failure under the configured memory budget'}
ref_200_label = state_labels.get(ref_200, ref_200)
ref_400_label = state_labels.get(ref_400, ref_400)
if (int(full_attempt.N) != 715512 or full_attempt.engine != 'SpatialCellChat' or
        full_attempt.reference_result_generated):
    raise SystemExit('Unexpected full-scale SpatialCellChat attempt record')
allocation = float(full_attempt.requested_allocation_GiB)
rss_gb = float(full_attempt.peak_RSS_GB_decimal)
wall_seconds = float(full_attempt.whole_process_seconds)

readme_path = ROOT/'README.md'
readme = readme_path.read_text()
readme = re.sub(
    r'(?:In D/E,|Panels D/E use).*?(?=\n## Main deliverables)',
    'Panels D/E use the original SpatialESS-blue / SpatialCellChat-orange encoding, with solid, '
    'dashed and dotted lines for CosMx, Atera breast and Atera cervical. Their shared legends '
    'identify method and cohort; paired-size shading is retained. Panel I reads the cervical 200k '
    'and 400k completion states from the validated extension run table. Its size ticks stand alone '
    'without an explanatory axis title. The 169,940-cell breast category is omitted from I only; '
    'its measurements remain in D/E and source data. Cervical SpatialESS completed both 200k and '
    '400k in the run table; the 400k category is now hidden in I, along with the nominal 450k '
    '(443,515-cell CosMx) category. Their run states remain in source data; the CosMx 443,515-cell '
    'measurements remain in D/E. The 700k reference attempt failed on a '
    f'{allocation:.1f} GiB dense intermediate request; the figure uses the same reference-failed '
    'cross as earlier failed attempts, and the exact failure is retained in source data. The visible 200k tick is '
    'followed directly by 700k.\n',
    readme, count=1, flags=re.S)
readme_path.write_text(readme)

qa_path = ROOT/'reports/QA_REPORT.md'
qa = qa_path.read_text()
qa = re.sub(
    r'Panels D and E encode methods.*?identical sizes\.',
    'Panels D and E use blue circles for SpatialESS and orange squares for SpatialCellChat. '
    'Solid, dashed and dotted lines distinguish CosMx, Atera breast and Atera cervical; shared '
    'legends identify both encodings. Neutral shading marks only input sizes with completed '
    'measurements from both methods.',
    qa, count=1, flags=re.S)
qa = re.sub(
    r'Panel I uses equally spaced .*?(?= The D, E and I panel previews)',
    f"Panel I uses equally spaced size ticks at 1k, 5k, 10k, 25k, 50k, 100k, 200k and 700k; "
    f"200k is followed directly by 700k and the explanatory x-axis title is omitted. The 400k "
    f"cervical and nominal 450k (443,515-cell CosMx) categories are hidden in I, with all run "
    f"states retained in source data; the 443,515-cell measurements remain in D/E. The "
    f"169,940-cell breast category and two 250-cell CosMx states are also hidden only in I. "
    f"At 700k, the same reference-failed cross marks the attempted run, which failed requesting "
    f"a {allocation:.1f} GiB dense intermediate under the 96 GiB virtual-memory limit; no "
    f"reference result was generated.",
    qa, count=1, flags=re.S)
qa = qa.replace(
    'The runtime axis includes all successful measurements, including the 2.39-6.02 s CosMx '
    'measurements below the former lower display limit.',
    'The restored D runtime axis starts at 8 s; smaller CosMx measurements remain in source '
    'data but fall below the displayed range.')
check_counts = {}
for name in ['data_checks', 'export_checks', 'source_data_checks']:
    checks = pd.read_csv(ROOT/'reports'/f'{name}.tsv', sep='\t')
    if not checks.status.eq('PASS').all():
        raise SystemExit(f'Failed validation in {name}')
    check_counts[name] = len(checks)
qa = re.sub(
    r'(?:This visual revision retained|The frozen baseline inputs and author workflow retain).*\Z',
    'The frozen baseline inputs and author workflow retain byte identity in the provenance '
    'checks; the cervical extension table was synchronized to terminal run records. '
    f"All {check_counts['data_checks']} numerical checks, "
    f"{check_counts['export_checks']} export checks and "
    f"{check_counts['source_data_checks']} baseline source checks passed. "
    'A color-vision-deficiency simulation was not performed.\n',
    qa, count=1, flags=re.S)
qa_path.write_text(qa)

print(f'QA report updated: SpatialESS 200k={ess_200}; 400k={ess_400}.')
