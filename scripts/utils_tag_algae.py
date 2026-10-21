"""Concatenate algal proteomes, prefixing each ID with its species (needed by 03_classify.py).
Reads data/raw/host_*.faa directly (stdin is ignored, kept for pipeline readability)."""
import glob, sys
from Bio import SeqIO
for f in sorted(glob.glob("data/raw/host_*.faa")):
    sp = f.split("host_")[1][:-4]
    for r in SeqIO.parse(f, "fasta"):
        sys.stdout.write(f">{sp}|{r.id}\n{r.seq}\n")
