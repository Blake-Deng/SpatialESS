"""Shared provenance, measurement and checkpoint helpers for the scale ladder."""
import csv
import hashlib
import json
import os
import re
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CFG = json.loads((ROOT / 'config/design.json').read_text())
BASE = Path(CFG['base_analysis'])
SIZES = CFG['paired_sizes'] + CFG['extension_sizes']


def stamp():
    return datetime.now().astimezone().isoformat()


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, default=str) + '\n')
    temporary.replace(path)


def sha256(path):
    value = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1048576), b''):
            value.update(block)
    return value.hexdigest()


def read_row(path):
    with path.open() as handle:
        return next(csv.DictReader(handle, delimiter='\t'))


def fixture(size):
    return 'full' if size == 'full' else 'central_' + str(size)


def tag(tissue, size, engine):
    return f'{tissue}_{fixture(size)}_{engine}_formal_b100'


def n_cells(tissue, size):
    if size != 'full':
        return int(size)
    return int(json.loads((BASE / 'inputs' / tissue / 'matrix_complete.json').read_text())['eligible_cells'])


def process_alive(pid):
    try:
        status = Path(f'/proc/{pid}/status').read_text()
    except FileNotFoundError:
        return False
    return not re.search(r'^State:\s+Z', status, re.M)


def time_metrics(path):
    if not path.exists():
        return {}
    text = path.read_text(errors='replace')
    result = {}
    match = re.search(r'Maximum resident set size \(kbytes\):\s+(\d+)', text)
    if match:
        result['peak_RSS_KiB'] = int(match[1])
        result['peak_RSS_GB_decimal'] = int(match[1]) * 1024 / 1e9
    match = re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s+([\d:.]+)', text)
    if match:
        parts = [float(x) for x in match[1].split(':')]
        result['whole_process_seconds'] = sum(v * 60 ** i for i, v in enumerate(reversed(parts)))
    match = re.search(r'Exit status:\s+(\d+)', text)
    if match:
        result['reported_exit_status'] = int(match[1])
    match = re.search(r'Command terminated by signal (\d+)', text)
    if match:
        result['termination_signal'] = int(match[1])
    return result


def failure_class(code, log_path, timed_out=False):
    if timed_out:
        return 'timeout'
    if code == 0:
        return 'completed'
    text = log_path.read_text(errors='replace') if log_path.exists() else ''
    if 'caught segfault' in text or code in (-11, 139):
        return 'segmentation_fault'
    if re.search(r'cannot allocate|std::bad_alloc|vector memory exhausted|Cannot allocate memory', text, re.I):
        return 'allocation_failure_under_memory_budget'
    if re.search(r'long vectors not supported|too many nonzero|32.bit|integer overflow', text, re.I):
        return 'implementation_limit_failure'
    if code in (-9, 137):
        return 'signal_9_unclassified'
    if code in (-6, 134):
        return 'abort'
    return 'implementation_error' if code is not None else 'missing_exit_evidence'


def archive_attempt(name, reason):
    destination = ROOT / 'logs/history' / (name + '_' + datetime.now().strftime('%Y%m%dT%H%M%S%f'))
    destination.mkdir(parents=True)
    for suffix in ['.exit.json', '.log', '.time.txt']:
        source = ROOT / 'logs' / (name + suffix)
        if source.exists():
            source.rename(destination / source.name)
    write_json(destination / 'archive_reason.json', {'reason': reason, 'time': stamp()})


def validate_success(name):
    result = ROOT / 'results' / name
    if not (result / 'complete.json').exists():
        return False
    row = read_row(result / 'summary.tsv')
    required = {'genes': 18028, 'lr': 3081, 'nboot': 100, 'seed_use': 1}
    if any(int(row[key]) != value for key, value in required.items()):
        raise RuntimeError('Saved result contract changed: ' + name)
    source = ROOT / 'inputs' / row['tissue'] / (row['fixture'] + '.rds')
    h = hashlib.md5()
    with source.open('rb') as handle:
        for block in iter(lambda: handle.read(1048576), b''):
            h.update(block)
    if h.hexdigest() != row['input_md5']:
        raise RuntimeError('Saved result input MD5 mismatch: ' + name)
    return True
