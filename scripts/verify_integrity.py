"""Verify the publication package against its recorded SHA-256 checksums."""
from pathlib import Path
import hashlib
root=Path(__file__).resolve().parents[1]
count=0
for line in (root/'SHA256SUMS').read_text().splitlines():
    expected,relative=line.split('  ',1)
    actual=hashlib.sha256((root/relative).read_bytes()).hexdigest()
    if actual!=expected: raise SystemExit(f'Checksum mismatch: {relative}')
    count+=1
print(f'PASS: {count} package files match SHA-256 checksums.')
