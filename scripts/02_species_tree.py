"""Global orthogroups (MMseqs2 clusters) -> single-copy near-universal markers ->
concatenated protein alignment -> IQ-TREE ML species tree."""
import subprocess, collections
from pathlib import Path
from Bio import SeqIO

OG = Path("work/og"); TREE = Path("work/tree"); TREE.mkdir(parents=True, exist_ok=True)
prot = {r.id: str(r.seq) for r in SeqIO.parse("work/genes/proteins.faa", "fasta")}
genomes = sorted({g.split("|")[0] for g in prot})
cl = collections.defaultdict(list)
for line in open(OG/"all_cluster.tsv"):
    rep, mem = line.split()
    cl[rep].append(mem)
# write orthogroup table (used downstream)
with open(OG/"orthogroups.tsv", "w") as fh:
    fh.write("og\tgene\n")
    for i, (rep, mems) in enumerate(sorted(cl.items(), key=lambda x: -len(x[1]))):
        for m in mems: fh.write(f"OG{i:05d}\t{m}\n")
ogs = collections.defaultdict(list)
for line in list(open(OG/"orthogroups.tsv"))[1:]:
    o, g = line.split(); ogs[o].append(g)
markers = []
for o, mems in ogs.items():
    c = collections.Counter(m.split("|")[0] for m in mems)
    if len(c) >= 0.95*len(genomes) and max(c.values()) == 1:
        markers.append(o)
print("genomes", len(genomes), "markers", len(markers))
concat = {g: [] for g in genomes}
for o in markers:
    fa = TREE/f"{o}.faa"
    SeqIO.write([SeqIO.SeqRecord.SeqRecord(SeqIO.Seq.Seq(prot[m]), id=m.split("|")[0], description="") for m in ogs[o]], fa, "fasta") if False else None
    with open(fa, "w") as fh:
        for m in ogs[o]: fh.write(f">{m.split('|')[0]}\n{prot[m]}\n")
    aln = subprocess.run(["mafft", "--auto", "--quiet", "--thread", "2", str(fa)], capture_output=True, text=True).stdout
    recs = {r.id: str(r.seq) for r in SeqIO.parse(__import__("io").StringIO(aln), "fasta")}
    L = len(next(iter(recs.values())))
    keep = [i for i in range(L) if sum(s[i] == "-" for s in recs.values()) <= 0.2*len(recs)]
    for g in genomes:
        s = recs.get(g, "-"*L)
        concat[g].append("".join(s[i] for i in keep))
with open(TREE/"concat.faa", "w") as fh:
    for g in genomes: fh.write(f">{g}\n{''.join(concat[g])}\n")
open(TREE/"markers.txt","w").write("\n".join(markers)+"\n")
