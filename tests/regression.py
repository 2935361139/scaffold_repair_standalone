#!/usr/bin/env python3
"""Real historical 3-sample/57-arm regression; no neural generation or checkpoint."""
import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def compare(fixture, output, exact_sdf=True):
    from scaffold_repair.evaluation import compare_endpoint
    from scaffold_repair.io import read_json, require
    old_root = fixture / 'intervention/group_01/rescue_v2'
    old = read_json(old_root / 'results.json')
    new = read_json(output / 'intervention/results.json')
    require(new['summary']['completed'] and not new['summary']['branch_errors'], 'New intervention run incomplete')
    expected = {(r['run_id'], r['arm']): r for r in old['records']}
    actual = {(r['run_id'], r['arm']): r for r in new['records']}
    wanted = {k for k in expected if k[1] != 'baseline'}
    require(wanted <= set(actual), 'Missing original intervention branches')
    sdf_equal = 0
    for key in sorted(wanted):
        a, b = actual[key], expected[key]
        compare_endpoint(a['evaluation'], b['evaluation'])
        require(a['receipt'] == b['receipt'], f'Optimizer/control receipt changed: {key}')
        rid, arm = key
        current = output / 'intervention/arms' / rid / (arm + '.sdf')
        previous = old_root / rid / (arm + '.sdf')
        equal = current.read_bytes() == previous.read_bytes()
        sdf_equal += equal
        if exact_sdf:
            require(equal, f'SDF bytes changed: {key}')
    golden_time = read_json(fixture / 'generation/group_01/localization_v2/geometry/temporal_evidence.json')
    times = {(r['id'], r['check'], r['stream'], r['sampler_atom_index'], r['environment_atom_index']): r for r in golden_time}
    time_count, pair_count = 0, 0
    for p in sorted((output / 'generated/samples').glob('*/diagnosis.json')):
        rid = p.parent.name
        diag = read_json(p)
        pair_count += len(diag['pairs'])
        for row in diag['temporal']['observations']:
            key = (rid, row['check'], row['stream'], row['sampler_atom_index'], row['environment_atom_index'])
            reference = times[key]
            for field in ('points', 'status', 'first_observed', 'first_persistent_observation', 'causal_onset'):
                require(row.get(field) == reference.get(field), f'Timing evidence changed: {key} {field}')
            time_count += 1
    require(time_count == len(times), 'Missing trajectory evidence')
    return dict(baseline_samples=3, intervention_control_endpoints_matched=len(wanted),
                optimizer_and_control_receipts_exact=len(wanted), sdf_files_byte_identical=sdf_equal,
                clash_pairs=pair_count, trajectory_series_exact=time_count,
                repaired=new['summary']['rescued'], repaired_with_relaxed_fixed_coordinates=new['summary']['rescued_with_relaxed_fixed_coordinates'],
                additional_v2_arms=len(actual)-len(wanted),
                note='Original 57 branch endpoints, receipts, SDF bytes and trajectory evidence must match; new v2 branches may change final selection.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--compare-only', action='store_true', help='Compare an already finished regression output')
    parser.add_argument('--fixture', type=Path, help='Use an already extracted reference bundle')
    args = parser.parse_args()
    out = args.output.resolve()
    if args.fixture:
        fixture = args.fixture.resolve()
    else:
        fixture = out.parent / (out.name + '_reference')
        if not args.compare_only:
            if fixture.exists(): raise ValueError('Reference extraction already exists; choose a new output name')
            fixture.mkdir(parents=True)
            with zipfile.ZipFile(ROOT / 'tests/data/regression3.zip') as archive:
                for entry in archive.infolist():
                    p = Path(entry.filename)
                    if p.is_absolute() or '..' in p.parts: raise ValueError('Unsafe archive entry')
                    archive.extract(entry, fixture)
    if not args.compare_only:
        cfg = out.parent / (out.name + '_config.py')
        if cfg.exists(): raise ValueError('Regression config already exists')
        cfg.parent.mkdir(parents=True, exist_ok=True)
        cfg.write_text('\n'.join(f'{key} = {value!r}' for key, value in dict(
            INPUT_MODE='flowr_records', FLOWR_GENERATED_DIR=str(fixture/'generation'),
            PROTEIN_PDB='', OUTPUT_DIR=str(out), INTERVENTION_INPUT_DIR=None, FAILURE_IDS=None).items()) + '\n')
        subprocess.run([sys.executable, str(ROOT/'run.py'), '--config', str(cfg)], check=True)
    result = compare(fixture, out)
    (out/'regression_comparison.json').write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
