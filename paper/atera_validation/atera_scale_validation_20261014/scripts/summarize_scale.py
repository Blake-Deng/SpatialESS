"""Export only measured values and distinguish completion from agreement."""
import json
import zipfile
import pandas as pd
from scale_common import BASE, CFG, ROOT, SIZES, fixture, n_cells, read_row, sha256, stamp, tag, time_metrics, write_json

STOP_STATES = {'timeout', 'allocation_failure_under_memory_budget', 'implementation_error',
               'segmentation_fault', 'implementation_limit_failure', 'abort', 'signal_9_unclassified',
               'missing_exit_evidence', 'failed_exit_1', 'failed_exit_3'}


def run_row(tissue, size, engine, baseline=False):
    name = tag(tissue, size, engine)
    source = BASE if baseline else ROOT
    result = source / 'results' / name
    reused = result.is_symlink()
    metrics_root = BASE if reused else source
    row = dict(tag=name, tissue=tissue, fixture=fixture(size), N=n_cells(tissue, size),
               engine=engine, phase='formal', genes=18028, lr=3081, nboot=100, seed_use=1,
               status='not_attempted', output_scope='historical_1k_5k_baseline' if baseline else 'extension',
               result_source=str(result.resolve()), reused_saved_result=reused or baseline,
               timeout_seconds=(3600 if baseline and engine == 'SpatialCellChat' else
                                CFG['reference_timeout_seconds'] if engine == 'SpatialCellChat' else CFG['ess_timeout_seconds']))
    exitfile = source / 'logs' / (name + '.exit.json')
    livefile = source / 'logs' / (name + '.running.json')
    if exitfile.exists():
        meta = json.loads(exitfile.read_text())
        row.update(status=meta['status'], returncode=meta['returncode'],
                   failure_log=str(source / 'logs' / (name + '.log')),
                   elapsed_until_exit_seconds=meta.get('elapsed_seconds'))
        for field in ['started', 'finished', 'adopted']:
            if field in meta:
                row[field] = meta[field]
        if 'observed_R_high_water_RSS_KiB' in meta:
            row['observed_R_high_water_RSS_KiB'] = meta['observed_R_high_water_RSS_KiB']
            row['observed_R_high_water_RSS_GB_decimal'] = meta['observed_R_high_water_RSS_KiB'] * 1024 / 1e9
        setup_log = source / 'logs' / (name + '.log')
        if (meta.get('returncode') == 1 and setup_log.exists() and
                'packageVersion("SpatialESS")' in setup_log.read_text()):
            row['status'] = 'setup_failure_excluded'
            row['failure_log'] = str(setup_log)
    if livefile.exists():
        live = json.loads(livefile.read_text())
        row.update(status='RUNNING', started=live['started'], heartbeat=live.get('heartbeat'),
                   observed_R_high_water_RSS_KiB=live.get('observed_R_high_water_RSS_KiB'))
    elif not baseline and (ROOT / 'reports/progress.json').exists():
        progress = json.loads((ROOT / 'reports/progress.json').read_text())
        if progress.get('tag') == name and progress.get('status') == 'RUNNING':
            row.update(status='RUNNING', started=progress['started'])
    if (result / 'complete.json').exists():
        summary = read_row(result / 'summary.tsv')
        for key in ['loading_seconds', 'inference_seconds', 'active_records', 'retained_records',
                    'graph_edges', 'input_md5', 'groups', 'raw_p_below_005']:
            row[key] = summary.get(key)
        if row['status'] in {'not_attempted', 'completed'}:
            row['status'] = 'completed'
    measures = time_metrics(metrics_root / 'logs' / (name + '.time.txt'))
    # Failed-attempt RSS/wall values stay identified as failed measurements.
    for field in ['peak_RSS_KiB', 'peak_RSS_GB_decimal', 'whole_process_seconds']:
        if field in measures:
            row[field if row['status'] == 'completed' else 'failed_attempt_' + field] = measures[field]
    if row['status'] in STOP_STATES:
        row['inference_seconds'] = None
        row['performance_note'] = 'Failed attempt; excluded from successful timing and speedup.'
    else:
        row['performance_note'] = 'One measurement; ESS reuse retains original measurement and hardware/threads.'
    return row


def summary():
    rows, comparisons = [], []
    for tissue in CFG['tissues']:
        for size in [1000, 5000] + SIZES:
            baseline = size in [1000, 5000]
            for engine in ['SpatialESS', 'SpatialCellChat']:
                rows.append(run_row(tissue, size, engine, baseline))
            path = ((BASE / 'results' / tag(tissue, size, 'SpatialESS')) if baseline else
                    (ROOT / 'comparisons' / tag(tissue, size, 'SpatialESS'))) / 'agreement.tsv'
            if path.exists():
                row = read_row(path)
                row['comparison_scope'] = 'historical_1k_5k_baseline' if baseline else 'new_extension'
                comparisons.append(row)
    perf = pd.DataFrame(rows)
    agreement = pd.DataFrame(comparisons)
    boundaries = []
    for tissue in CFG['tissues']:
        ref = perf[(perf.tissue == tissue) & (perf.engine == 'SpatialCellChat')].sort_values('N')
        success = ref[ref.status == 'completed']
        failed = ref[ref.status.isin(STOP_STATES)]
        passed = agreement[(agreement.tissue == tissue) & agreement.strict_pass.astype(str).str.lower().eq('true')]
        last = int(success.N.max()) if len(success) else None
        first = int(failed.iloc[0].N) if len(failed) else None
        first_state = failed.iloc[0].status if len(failed) else None
        fulln = n_cells(tissue, 'full')
        at_full = last == fulln
        alltested = ref.status.eq('completed').all()
        boundaries.append(dict(
            tissue=tissue, maximum_reference_completed_N=last,
            maximum_strictly_validated_N=int(pd.to_numeric(passed.N).max()) if len(passed) else None,
            first_reference_uncompleted_N=first, first_uncompleted_state=first_state,
            current_running_N=int(ref[ref.status == 'RUNNING'].N.iloc[0]) if ref.status.eq('RUNNING').any() else None,
            remaining_not_attempted_scales=','.join(ref[ref.status.isin(['not_attempted', 'setup_failure_excluded'])].fixture),
            full_tissue_reference_completed=at_full, actual_performance_limit_established=False,
            interpretation=('All predefined scales completed; no tested upper limit found.' if alltested else
                            'Completion boundary under the fixed 12-hour/96-GiB virtual-address budget; not a universal engine limit.' if first else
                            'Ladder in progress; no failure boundary established.')))
    boundaries = pd.DataFrame(boundaries)
    paired = []
    for row in comparisons:
        tissue, name = row['tissue'], row['fixture']
        selected = perf[(perf.tissue == tissue) & (perf.fixture == name) & (perf.status == 'completed')]
        if len(selected) != 2:
            continue
        ess, ref = selected[selected.engine == 'SpatialESS'].iloc[0], selected[selected.engine == 'SpatialCellChat'].iloc[0]
        if not pd.isna(ess.get('whole_process_seconds')) and not pd.isna(ref.get('whole_process_seconds')):
            paired.append(dict(tissue=tissue, fixture=name, N=int(row['N']),
                               whole_process_speedup=float(ref.whole_process_seconds) / float(ess.whole_process_seconds),
                               reference_over_ESS_peak_RSS=float(ref.peak_RSS_GB_decimal) / float(ess.peak_RSS_GB_decimal),
                               strict_pass=str(row['strict_pass']).lower() == 'true',
                               interpretation='Same prepared fixture; reused ESS timing is explicitly marked.'))
    tables = ROOT / 'tables'
    tables.mkdir(exist_ok=True)
    frames = {'atera_scale_performance': perf, 'atera_scale_agreement': agreement,
              'atera_scale_boundaries': boundaries, 'atera_scale_paired_ratios': pd.DataFrame(paired)}
    for name, frame in frames.items():
        for delimiter, ext in [('\t', 'tsv'), (',', 'csv')]:
            frame.to_csv(tables / (name + '.' + ext), sep=delimiter, index=False)
    temporary = tables / 'Atera_scale_validation.tmp.xlsx'
    with pd.ExcelWriter(temporary, engine='openpyxl') as writer:
        for name, frame in frames.items():
            frame.to_excel(writer, sheet_name=name.replace('atera_scale_', ''), index=False)
            writer.sheets[name.replace('atera_scale_', '')].freeze_panes = 'A2'
        pd.DataFrame([
            ['Input', '18,028 measured genes; 3,081 fixed LR; 100 permutations; seed=1. No change to filters.'],
            ['Scope', 'Fixed-LR benchmark, not standard overexpression LR-selection workflow.'],
            ['CPU', f'Logical CPU {CFG["cpu"]}; single-thread BLAS/OpenMP; serial workers.'],
            ['Budget', 'Reference 12 h per attempt; 96 GiB virtual-address limit, not measured RSS.'],
            ['Failure', 'Timeout/allocation/segfault/error saved separately; signal 9 is not called OOM without evidence.'],
            ['Performance', 'Single measurements. Successful reference times only; failed times not used for speedups.'],
            ['Completion', 'Maximum completed and first uncompleted under budget; no claim of intrinsic performance limit.'],
            ['Agreement', 'Full keys, probability max error <=1e-15, exact p, both <.05 and <=.05 decisions.'],
            ['Reuse', '25k/50k/100k/full ESS results read only; new comparisons written independently.'],
            ['Biology', 'No new inference definition, FDR, biological discovery or causal validation.'],
        ], columns=['Item', 'Definition']).to_excel(writer, sheet_name='definitions', index=False)
    temporary.replace(tables / 'Atera_scale_validation.xlsx')
    text = ['# Atera reference scale validation', '', 'Updated: ' + stamp(), '',
            'The table distinguishes successful execution, numerical agreement, and budget-limited failures.', '',
            '| Tissue | Maximum reference completed N | Maximum strict N | First uncompleted N | State |',
            '| --- | ---: | ---: | ---: | --- |']
    for row in boundaries.to_dict('records'):
        text.append('| ' + ' | '.join('NA' if pd.isna(row[k]) else str(row[k]) for k in [
            'tissue', 'maximum_reference_completed_N', 'maximum_strictly_validated_N',
            'first_reference_uncompleted_N', 'first_uncompleted_state']) + ' |')
    text += ['', 'Unattempted or running scales do not establish an engine limit. Each tissue stops at its first failed reference run or agreement discrepancy. The other tissue continues its predefined ladder.',
             '', 'Old 1k/5k results used a one-hour reference timeout; extension attempts use twelve hours. CPU, thread count, virtual address cap and scientific inputs are held fixed. ESS timings reused from the original experiment retain their source.',
             '', 'Raw permutation p < 0.05 is descriptive; no new FDR correction or biological validation is added.']
    (ROOT / 'reports/Atera_scale_results_summary.md').write_text('\n'.join(text) + '\n')
    write_json(ROOT / 'reports/table_status.json', dict(updated=stamp(),
        completed_reference_extension=int(((perf.engine == 'SpatialCellChat') & (perf.output_scope == 'extension') & (perf.status == 'completed')).sum()),
        strict_extension_pairs=int((agreement.comparison_scope.eq('new_extension') & agreement.strict_pass.astype(str).str.lower().eq('true')).sum()),
        boundaries=boundaries.astype(object).where(boundaries.notna(), None).to_dict('records')))
    progress = ROOT / 'reports/progress.json'
    if progress.exists() and json.loads(progress.read_text()).get('status') in {
            'COMPLETE_TESTED_PLAN', 'FINISHED_WITH_BUDGET_OR_VALIDATION_STOP'}:
        archive = ROOT / 'Atera_reference_scale_validation_results.zip'
        with zipfile.ZipFile(archive, 'w', zipfile.ZIP_DEFLATED) as handle:
            for folder in ['tables', 'reports', 'logs', 'config', 'scripts', 'comparisons']:
                for path in sorted((ROOT / folder).rglob('*')):
                    if path.is_file() and '__pycache__' not in path.parts:
                        handle.write(path, str(path.relative_to(ROOT)))
            if (ROOT / 'README.md').exists():
                handle.write(ROOT / 'README.md', 'README.md')
            for path in (ROOT / 'results').glob('*/summary.tsv'):
                handle.write(path, str(path.relative_to(ROOT)))
        with zipfile.ZipFile(archive) as handle:
            if handle.testzip() is not None:
                raise RuntimeError('Final ZIP CRC failed')
        write_json(ROOT / 'reports/archive_verification.json', dict(
            archive=str(archive), bytes=archive.stat().st_size, sha256=sha256(archive), crc='PASS'))


if __name__ == '__main__':
    summary()
