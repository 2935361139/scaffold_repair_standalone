#!/usr/bin/env python3
"""Single configurable entrance for external molecular data; no repair or labels."""
import argparse
import runpy
from pathlib import Path
from universal_input.normalize import convert

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=Path, default=Path(__file__).with_name('import_config.py'))
    args = parser.parse_args()
    path = args.config.resolve()
    result = convert(runpy.run_path(str(path))['IMPORT'], path.parent)
    print('IMPORTED:', result['total_samples'], 'samples')
    print('STATUS:', result['status'])
    print('REPORT:', result['report'])
    for item in result['needs']:
        print('NEEDED:', item)
    import zipfile
    output = Path(result['report']).parent
    package = output / 'conversion_review.zip'
    files = [p for p in output.rglob('*') if p.is_file() and p != package]
    with zipfile.ZipFile(package, 'x', zipfile.ZIP_DEFLATED) as z:
        for p in sorted(files):
            z.write(p, p.relative_to(output))
    print('SEND THIS FILE:', package)
