"""Refresh distributed-file hashes after deliberate documentation changes.
This does not modify SOURCE_PROVENANCE.json or any scientific file.
"""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parent
EXCLUDED={'.git','.venv','runs','__pycache__'}
manifest={}
for p in sorted(ROOT.rglob('*')):
    relative=p.relative_to(ROOT)
    if not p.is_file() or any(part in EXCLUDED for part in relative.parts):continue
    if p.name=='MANIFEST_SHA256.json' or p.suffix in ['.pyc','.pyo']:continue
    manifest[relative.as_posix()]=hashlib.sha256(p.read_bytes()).hexdigest()
(ROOT/'MANIFEST_SHA256.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
print('Manifest updated for',len(manifest),'files.')
