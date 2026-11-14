#!/usr/bin/env python3
"""Check already-generated molecules, then run the unchanged finite repair bank."""
import argparse
import importlib.util
import os
from pathlib import Path
import sys

# Match the single-thread numerical execution used in the recorded experiments.
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=Path(__file__).with_name('config.py'))
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--check-only', '--generate-only', dest='check_only', action='store_true',
                       help='Import/check/classify existing molecules; never call a generator')
    modes.add_argument('--intervene-only', action='store_true', help='Read checked failures; do not regenerate or reselect inputs')
    args = parser.parse_args(argv)
    from scaffold_repair.io import require, resolve
    from scaffold_repair.pipeline import check, intervene
    path = args.config.resolve()
    spec = importlib.util.spec_from_file_location('user_experiment_config', path)
    config = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(config)
    output = resolve(config.OUTPUT_DIR, path.parent)
    if not args.intervene_only:
        if config.INPUT_MODE == 'flowr_records':
            from scaffold_repair.adapters.flowr_records import load
            protein = resolve(config.PROTEIN_PDB, path.parent) if config.PROTEIN_PDB else None
            items = load(resolve(config.FLOWR_GENERATED_DIR, path.parent), protein)
        elif config.INPUT_MODE == 'sdf':
            from scaffold_repair.adapters.sdf import load_sdf
            require(config.PROTEIN_PDB, 'Fill PROTEIN_PDB for SDF input')
            items = load_sdf(*(resolve(getattr(config, name), path.parent) for name in
                ('INPUT_SDF', 'PROTEIN_PDB', 'REFERENCE_SDF', 'ATOM_MAPPING_FILE')),
                record_indices=config.INPUT_RECORD_INDICES)
        elif config.INPUT_MODE == 'manifest':
            from scaffold_repair.adapters.sdf import load_manifest
            items = load_manifest(resolve(config.MANIFEST_PATH, path.parent))
        else:
            raise ValueError('INPUT_MODE must be flowr_records, sdf or manifest')
        result = check(items, output / 'generated')
        print(f"CHECK COMPLETE: {result['evaluated_samples']} total; {result['passed']} passed; {result['failed']} failed")
        print(f"RESULT: {output / 'generated/review.zip'}")
    if not args.check_only:
        failed = resolve(config.INTERVENTION_INPUT_DIR, path.parent) if config.INTERVENTION_INPUT_DIR else output / 'generated/failed'
        result = intervene(failed, output / 'intervention', config.FAILURE_IDS)
        print(f"INTERVENTION COMPLETE: {result['rescued']}/{result['original_failed']} original failures repaired; "
              f"{result['rescued_with_relaxed_fixed_coordinates']} required relaxed fixed coordinates; "
              f"{result['branch_errors']} unevaluated branches")
        print(f"RESULT: {output / 'intervention/review.zip'}")


if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        print(f'STOPPED: {type(exc).__name__}: {exc}', file=sys.stderr)
        sys.exit(1)
