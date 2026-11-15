"""Install two small read/write hooks; back up, validate version, never overwrite unknown edits."""
from pathlib import Path
import argparse
import hashlib
from datetime import datetime
import shutil


def install(repo):
    changes={}
    rules={
        'scaffold_repair/evaluation.py': [
            ('from .io import require', 'from .io import require\nfrom .import_state import restore as restore_import_state, encode as encode_import_state, PROPERTY as IMPORT_STATE_PROPERTY'),
            ('    return supplier[0]', '    return restore_import_state(supplier[0])'),
            ('            writer.write(mol)', '            if mol.HasProp(IMPORT_STATE_PROPERTY):\n                writer.SetKekulize(False)\n                writer.SetForceV3000(True)\n                encode_import_state(mol)\n            writer.write(mol)')],
        'scaffold_repair/adapters/sdf.py': [
            ('from ..evaluation import task_from_mapping', 'from ..evaluation import task_from_mapping\nfrom ..import_state import restore as restore_import_state'),
            ('    mol = supplier[index]', '    mol = restore_import_state(supplier[index])'),
            ('    ref = refs[reference_index]', '    ref = restore_import_state(refs[reference_index])')]
    }
    expected = {'scaffold_repair/evaluation.py': 'd5de61e33474a66635ee668a9956390b364f39837f217c5667767a878c03f323', 'scaffold_repair/adapters/sdf.py': '22a50c4ae2dc8069b900dfb5091a017e93356309f71195525b15b6465cb2e2f5'}
    for rel, replacements in rules.items():
        p=repo/rel;text=p.read_text()
        if all(new in text for old,new in replacements):continue
        if hashlib.sha256(p.read_bytes()).hexdigest()!=expected[rel]:
            raise RuntimeError('Different local version; no changes made. Send: '+str(p))
        for old,new in replacements:
            if text.count(old)!=1:raise RuntimeError('Unexpected source layout: '+rel)
            text=text.replace(old,new)
        compile(text,str(p),'exec');changes[rel]=text
    if not changes:
        print('Already installed.');return
    backup=repo/'archive'/('input_adapter_'+datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    for rel in changes:
        target=backup/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(repo/rel,target)
    for rel,text in changes.items():(repo/rel).write_text(text)
    print('Installed. Backup:',backup)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--repo',type=Path,default=Path(__file__).resolve().parent)
    install(parser.parse_args().repo.resolve())
