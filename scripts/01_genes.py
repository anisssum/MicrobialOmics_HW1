"""Uniform gene calling (Prodigal, single mode, code 11) for all chlorovirus genomes.
Writes work/genes/{proteins.faa,cds.fna,upstream.fna,genes.tsv}."""
import subprocess, sys, os, re
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

RAW = Path("data/raw/chlorovirus_genomes.fna")
OUT = Path("work/genes"); OUT.mkdir(parents=True, exist_ok=True)
UP = 150  # nt upstream of start codon kept for promoter analysis
tmp = OUT / "tmp"; tmp.mkdir(exist_ok=True)

faa, fna, upf, rows = [], [], [], []
for rec in SeqIO.parse(RAW, "fasta"):
    acc = rec.id
    seq = str(rec.seq).upper()
    nfrac = seq.count("N") / len(seq)
    gpath = tmp / f"{acc}.fna"
    SeqIO.write(rec, gpath, "fasta")
    gff = tmp / f"{acc}.gff"
    subprocess.run(["prodigal", "-i", str(gpath), "-f", "gff", "-o", str(gff), "-g", "11", "-q"], check=True)
    n = 0
    for line in open(gff):
        if line.startswith("#"): continue
        f = line.rstrip().split("\t")
        s, e, strand = int(f[3]), int(f[4]), f[6]
        attrs = dict(kv.split("=") for kv in f[8].strip(";").split(";"))
        partial = attrs.get("partial", "00")
        if partial != "00": continue
        n += 1
        gid = f"{acc}|g{n:04d}"
        nt = Seq(seq[s-1:e])
        if strand == "-":
            nt = nt.reverse_complement()
            up = Seq(seq[e:e+UP]).reverse_complement()
        else:
            up = Seq(seq[max(0, s-1-UP):s-1])
        aa = nt.translate(table=11, to_stop=True)
        if len(aa) < 50: continue  # skip tiny ORFs (<50 aa)
        fna.append(f">{gid}\n{nt}\n"); faa.append(f">{gid}\n{aa}\n"); upf.append(f">{gid}\n{up}\n")
        rows.append(f"{gid}\t{acc}\t{s}\t{e}\t{strand}\t{len(aa)}")
    print(acc, len(seq), f"N={nfrac:.4f}", n, file=sys.stderr)

open(OUT/"proteins.faa","w").write("".join(faa))
open(OUT/"cds.fna","w").write("".join(fna))
open(OUT/"upstream.fna","w").write("".join(upf))
open(OUT/"genes.tsv","w").write("gene\tgenome\tstart\tend\tstrand\taa_len\n"+"\n".join(rows)+"\n")
