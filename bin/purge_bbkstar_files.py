from pathlib import Path

p = Path().absolute();

SUFFIXES_TO_KEEP = {".py", ".txt", ".stdout"}

for g in p.glob('*'):
    if g.is_file():
        if g.suffix not in SUFFIXES_TO_KEEP:
            g.unlink();