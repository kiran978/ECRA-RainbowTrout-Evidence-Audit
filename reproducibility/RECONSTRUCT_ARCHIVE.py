"""Restore private ECRA data snapshot from the GitHub part files.

Requires Python standard library only. Run from repository root: 
python reproducibility/RECONSTRUCT_ARCHIVE.py
"""
from pathlib import Path
import hashlib

root = Path(__file__).resolve().parent
parts = sorted((root/'parts').glob('ECRA_Complete_Private_Deposit_20261122.zip.part-*'))
if len(parts) != 15:
    raise SystemExit(f'Expected 15 parts; found {len(parts)}')
output = root/'ECRA_Complete_Private_Deposit_20261122.zip'
with output.open('wb') as dst:
    for part in parts:
        with part.open('rb') as f:
            while block := f.read(1024*1024):
                dst.write(block)
h = hashlib.sha256(output.read_bytes()).hexdigest()
expected='4e1836c317d3f882f374d53077b84291b109e2ff049895e000952d0d5c6ce080'
if h != expected:
    output.unlink(missing_ok=True)
    raise SystemExit(f'SHA-256 MISMATCH: expected {expected}; got {h}')
print(f'Archive rebuilt and SHA-256 verified: {output} ({output.stat().st_size} bytes)')
