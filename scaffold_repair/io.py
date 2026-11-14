"""Small, explicit on-disk records; no model code or saved Python objects."""
import csv
import hashlib
import json
import re
from pathlib import Path


def require(ok, message):
    if not ok:
        raise ValueError(message)


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + '.tmp')
    temp.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
    temp.replace(path)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def safe_id(value):
    require(isinstance(value, str) and re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,179}', value),
            f'Invalid sample ID: {value!r}. Use letters, digits, underscore, dot or hyphen.')
    return value


def resolve(path, base):
    return (Path(base) / Path(path).expanduser()).resolve()


def save_csv(path, rows):
    fields = list(dict.fromkeys(key for row in rows for key in row)) or ['status']
    with Path(path).open('w', newline='', encoding='utf-8-sig') as handle:
        writer = csv.DictWriter(handle, fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: json.dumps(value, ensure_ascii=False) if isinstance(value, (list, dict))
                             else value for key, value in row.items()})
