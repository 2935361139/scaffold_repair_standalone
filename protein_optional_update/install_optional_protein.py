"""Validated, backed-up installation. Does not touch config, inputs or experiment outputs."""
import argparse
from datetime import datetime
import hashlib
import json
from pathlib import Path
import shutil


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def inside(base, relative):
    path = (base / relative).resolve()
    if not path.is_relative_to(base):
        raise ValueError('Path outside project: ' + relative)
    return path


def install(repo):
    root = Path(__file__).resolve().parent
    manifest = json.loads((root/'package_manifest.json').read_text())
    if not (repo/'convert_inputs.py').is_file() or not (repo/'import_config.py').is_file():
        raise RuntimeError('Install the previous universal_input_update first. No files changed.')
    changes = []
    errors = []
    for rel, expected in manifest['unchanged_dependencies'].items():
        path = inside(repo, rel)
        if not path.is_file() or sha(path) != expected:
            errors.append('Different or missing dependency: ' + rel)
    for rel, versions in manifest['files'].items():
        source, target = inside(root/'payload', rel), inside(repo, rel)
        if sha(source) != versions['after']:
            errors.append('Damaged update package: ' + rel)
            continue
        compile(source.read_text(encoding='utf-8'), str(source), 'exec')
        actual = sha(target) if target.is_file() else None
        if actual == versions['after']:
            continue
        if actual != versions['before'] or (target.exists() and not target.is_file()):
            errors.append('Different local version: ' + rel)
        else:
            changes.append((rel, source, target))
    if errors:
        raise RuntimeError('No files changed. Preserve your files and send this output:\n' + '\n'.join(errors))
    if not changes:
        print('Already installed; no files changed.')
        return
    backup = repo/'archive'/('protein_optional_' + datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    backup.mkdir(parents=True)
    existed = {}
    for rel, source, target in changes:
        existed[rel] = target.is_file()
        if existed[rel]:
            dest = backup/rel; dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(target, dest)
    (backup/'installation.json').write_text(json.dumps(dict(files=existed, package_manifest=manifest), indent=2))
    written = []
    try:
        for rel, source, target in changes:
            target.parent.mkdir(parents=True, exist_ok=True)
            staged = target.with_name(target.name + '.protein_optional.tmp')
            staged.write_bytes(source.read_bytes())
            staged.replace(target)
            written.append((rel, target))
    except Exception:
        for rel, target in reversed(written):
            if existed[rel]: shutil.copy2(backup/rel, target)
            else: target.unlink()
        raise
    print('Installed:', len(changes), 'files')
    print('BACKUP:', backup)
    print('Your config.py, import_config.py, inputs and outputs were preserved.')
    print('Next: edit import_config.py as described in README_CN.md; use new input/output directories.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, required=True)
    try:
        install(parser.parse_args().repo.resolve())
    except Exception as error:
        raise SystemExit(str(error))
