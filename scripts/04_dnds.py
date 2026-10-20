"""Per-orthogroup dN/dS (PAML codeml, M0 one-ratio, F3x4) within each host-range clade.

For every (clade, OG) with >= MIN_SEQ single-copy, full-length members:
  protein MAFFT alignment -> codon back-translation -> IQ-TREE (GTR+G4, -fast) tree -> codeml M0.
Members are subsampled to <= MAX_SEQ sequences (fixed seed) to bound run time.
A second codeml run with omega fixed to 1 gives the neutral-evolution LRT.
Output: work/dnds/dnds_raw.tsv (resumable) and codon alignments in work/dnds/aln/.
"""
import os, re, sys, io, random, subprocess, tempfile, shutil
from multiprocessing import Pool
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

MIN_SEQ, MAX_SEQ, SEED = 6, 10, 42
OUT = Path("work/dnds"); (OUT / "aln").mkdir(parents=True, exist_ok=True)
RAW = OUT / "dnds_raw.tsv"
CLADES = sys.argv[1].split(",") if len(sys.argv) > 1 else ["NC64A", "SAG", "Pbi", "Osy"]

genes = pd.read_csv("work/genes_annotated.tsv", sep="\t")
cds = {r.id: str(r.seq).upper() for r in SeqIO.parse("work/genes/cds.fna", "fasta")}

CTL = """seqfile = codon.fa
treefile = tree.nwk
outfile = mlc
noisy = 0
verbose = 0
runmode = 0
seqtype = 1
CodonFreq = 2
clock = 0
model = 0
NSsites = 0
icode = 0
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = 0.4
cleandata = 1
"""


def build_jobs():
    jobs = []
    for (clade, og), d in genes.groupby(["clade", "og"]):
        if clade not in CLADES:
            continue
        cnt = d.genome.value_counts()
        d = d[d.genome.isin(cnt[cnt == 1].index)]          # single-copy genomes only
        med = d.aa_len.median()
        d = d[(d.aa_len >= 0.8 * med) & (d.aa_len <= 1.2 * med)]  # drop truncated / fused
        ids = [g for g in d.gene if cds[g].count("N") == 0]
        if len(ids) < MIN_SEQ:
            continue
        random.Random(f"{SEED}{og}{clade}").shuffle(ids)
        jobs.append((clade, og, sorted(ids[:MAX_SEQ])))
    return jobs


def codon_align(ids, wd):
    prot = "".join(f">s{i}\n{Seq(cds[g][:len(cds[g]) // 3 * 3]).translate(table=11).rstrip('*')}\n" for i, g in enumerate(ids))
    (wd / "p.fa").write_text(prot)
    aln = subprocess.run(["mafft", "--auto", "--quiet", "--thread", "1", "p.fa"], cwd=wd,
                         capture_output=True, text=True).stdout
    out = {}
    for r in SeqIO.parse(io.StringIO(aln), "fasta"):
        i = int(r.id[1:]); nt = cds[ids[i]]; k = 0; s = []
        for a in str(r.seq):
            if a == "-":
                s.append("---")
            else:
                s.append(nt[3 * k:3 * k + 3]); k += 1
        out[r.id] = "".join(s)
    return out


def run(job):
    clade, og, ids = job
    wd = Path(tempfile.mkdtemp(prefix=f"{og}_{clade}_", dir="/tmp"))
    try:
        aln = codon_align(ids, wd)
        L = len(next(iter(aln.values())))
        # remove codon columns containing stop codons (cleandata would drop them anyway)
        with open(wd / "codon.fa", "w") as fh:
            fh.write(f" {len(aln)} {L}\n")
            for k, v in aln.items():
                fh.write(f"{k}\n{v}\n")
        with open(OUT / "aln" / f"{clade}_{og}.fna", "w") as fh:
            for k, v in aln.items():
                fh.write(f">{ids[int(k[1:])]}\n{v}\n")
        with open(wd / "codon_iq.fa", "w") as fh:
            for k, v in aln.items():
                fh.write(f">{k}\n{v}\n")
        subprocess.run(["iqtree", "-s", "codon_iq.fa", "-m", "GTR+G4", "-fast", "-T", "1", "--quiet",
                        "-seed", "42", "--prefix", "t"], cwd=wd, capture_output=True)
        tree = (wd / "t.treefile").read_text().strip()
        tree = re.sub(r":[0-9.eE-]+", "", tree)  # topology only
        (wd / "tree.nwk").write_text(f" {len(aln)} 1\n{tree}\n")
        (wd / "codeml.ctl").write_text(CTL)
        subprocess.run(["codeml", "codeml.ctl"], cwd=wd, capture_output=True, timeout=900)
        mlc = (wd / "mlc").read_text()
        g = lambda pat: float(re.search(pat, mlc).group(1))
        # null model: neutral evolution, omega fixed to 1 (LRT, 1 df)
        (wd / "codeml.ctl").write_text(CTL.replace("fix_omega = 0", "fix_omega = 1").replace("omega = 0.4", "omega = 1").replace("outfile = mlc", "outfile = mlc1"))
        subprocess.run(["codeml", "codeml.ctl"], cwd=wd, capture_output=True, timeout=900)
        lnL1 = float(re.search(r"lnL\(.*?\):\s*(-?[0-9.]+)", (wd / "mlc1").read_text()).group(1))
        res = dict(clade=clade, og=og, n_seq=len(ids), n_codons=L // 3,
                   omega=g(r"omega \(dN/dS\)\s*=\s*([0-9.]+)"),
                   kappa=g(r"kappa \(ts/tv\)\s*=\s*([0-9.]+)"),
                   tree_dN=g(r"tree length for dN:\s*([0-9.]+)"),
                   tree_dS=g(r"tree length for dS:\s*([0-9.]+)"),
                   lnL=g(r"lnL\(.*?\):\s*(-?[0-9.]+)"), lnL_w1=lnL1)
    except Exception as e:
        res = dict(clade=clade, og=og, n_seq=len(ids), error=str(e)[:100])
    finally:
        shutil.rmtree(wd, ignore_errors=True)
    return res


if __name__ == "__main__":
    jobs = build_jobs()
    done = set()
    if RAW.exists():
        prev = pd.read_csv(RAW, sep="\t")
        done = set(zip(prev.clade, prev.og))
    jobs = [j for j in jobs if (j[0], j[1]) not in done]
    print(f"{len(jobs)} jobs to run", file=sys.stderr)
    header = not RAW.exists()
    with Pool(int(os.environ.get("THREADS", 2))) as p, open(RAW, "a") as fh:
        cols = ["clade", "og", "n_seq", "n_codons", "omega", "kappa", "tree_dN", "tree_dS", "lnL", "lnL_w1", "error"]
        if header:
            fh.write("\t".join(cols) + "\n")
        for i, r in enumerate(p.imap_unordered(run, jobs)):
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n"); fh.flush()
            if i % 50 == 0:
                print(i, file=sys.stderr)
