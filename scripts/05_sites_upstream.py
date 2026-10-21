"""Within-gene conservation landscape and upstream (promoter) constraint.

For each (clade, OG) codon alignment produced by 04_dnds.py:
 1. per-site amino-acid diversity (pi_aa, mean pairwise mismatch over non-gap pairs);
 2. pairwise synonymous divergence dS (Nei-Gojobori 1986) between members;
 3. heterogeneity of conservation along the protein: variance of 15-codon window means
    vs. 999 column permutations  -> p_hetero  ("are there conserved blocks?");
 4. for algal-like OGs: pi_aa inside vs. outside the region homologous to the algal protein
    (DIAMOND query coordinates mapped to alignment columns);
 5. upstream intergenic region (<=150 nt, truncated at the neighbouring gene): pairwise
    JC distance d_up, constraint index = median(d_up / dS) over pairs with 0.02 < dS < 1.5.
Output: work/sites/site_stats.tsv
"""
import io, os, sys, math, random, subprocess, tempfile
from multiprocessing import Pool
from pathlib import Path
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds

WIN, NPERM, MAXUP = 15, 999, 150
ALN = Path("work/dnds/aln"); OUT = Path("work/sites"); OUT.mkdir(parents=True, exist_ok=True)

genes = pd.read_csv("work/genes_annotated.tsv", sep="\t").set_index("gene")
genomes = {r.id: str(r.seq).upper() for r in SeqIO.parse("data/raw/chlorovirus_genomes.fna", "fasta")}

# upstream intergenic coordinates
up_seq = {}
for gnm, d in genes.reset_index().groupby("genome"):
    d = d.sort_values("start")
    starts, ends = d.start.values, d.end.values
    for gid, s, e, strand in zip(d.gene, d.start, d.end, d.strand):
        seq = genomes[gnm]
        if strand == "+":
            prev_end = ends[ends < s].max() if (ends < s).any() else 0
            a = max(prev_end, s - 1 - MAXUP); up_seq[gid] = seq[a:s - 1]
        else:
            nxt = starts[starts > e].min() if (starts > e).any() else len(seq) + 1
            b = min(nxt - 1, e + MAXUP); up_seq[gid] = str(Seq(seq[e:b]).reverse_complement())


def pi_aa_columns(prots):
    L = len(prots[0]); n = len(prots)
    pis, gapf = np.full(L, np.nan), np.zeros(L)
    for i in range(L):
        col = [p[i] for p in prots]
        aa = [c for c in col if c not in "-X*"]
        gapf[i] = 1 - len(aa) / n
        if len(aa) < 2:
            continue
        m = len(aa); same = sum(v * (v - 1) for v in pd.Series(aa).value_counts().values)
        pis[i] = 1 - same / (m * (m - 1))
    return pis, gapf


def jc(a, b):
    pairs = [(x, y) for x, y in zip(a, b) if x in "ACGT" and y in "ACGT"]
    if len(pairs) < 30:
        return np.nan
    p = sum(x != y for x, y in pairs) / len(pairs)
    return np.nan if p >= 0.74 else -0.75 * math.log(1 - 4 * p / 3)


def mafft(seqs):
    fa = "".join(f">s{i}\n{s}\n" for i, s in enumerate(seqs))
    r = subprocess.run(["mafft", "--auto", "--quiet", "--thread", "1", "-"], input=fa, capture_output=True, text=True).stdout
    return [str(x.seq).upper() for x in sorted(SeqIO.parse(io.StringIO(r), "fasta"), key=lambda x: int(x.id[1:]))]


def analyse(path):
    clade, og = path.stem.split("_", 1)
    recs = list(SeqIO.parse(path, "fasta"))
    ids = [r.id for r in recs]; nts = [str(r.seq) for r in recs]
    prots = ["".join("-" if s[i:i + 3] == "---" else str(Seq(s[i:i + 3]).translate(table=11)) for i in range(0, len(s), 3)) for s in nts]
    pis, gapf = pi_aa_columns(prots)
    ok = (gapf <= 0.2) & ~np.isnan(pis)
    res = dict(clade=clade, og=og, n_seq=len(ids), n_cols=int(ok.sum()), pi_aa=float(np.nanmean(pis[ok])) if ok.any() else np.nan)
    # heterogeneity of conservation
    v = pis[ok]
    if len(v) >= 3 * WIN and v.std() > 0:
        def winvar(x):
            c = np.convolve(x, np.ones(WIN) / WIN, mode="valid"); return c.var()
        obs = winvar(v); rng = np.random.default_rng(42)
        perm = np.array([winvar(rng.permutation(v)) for _ in range(NPERM)])
        res.update(hetero_ratio=obs / perm.mean(), p_hetero=(1 + (perm >= obs).sum()) / (NPERM + 1))
        # most conserved window (relative position along the alignment)
        c = np.convolve(v, np.ones(WIN) / WIN, mode="valid")
        res.update(min_window_pi=float(c.min()), max_window_pi=float(c.max()))
    # homologous vs non-homologous region (algal-like genes)
    sub = genes.loc[ids]
    if sub.algal_hit.mean() >= 0.5:
        cover = np.zeros(len(prots[0]))
        for p, (gid, row) in zip(prots, sub.iterrows()):
            if not row.algal_hit:
                continue
            k = 0
            for i, a in enumerate(p):
                if a != "-":
                    k += 1
                    if row.qstart <= k <= row.qend:
                        cover[i] += 1
        hom = cover >= 0.5 * sub.algal_hit.sum()
        a, b = pis[ok & hom], pis[ok & ~hom]
        if len(a) >= 20 and len(b) >= 20:
            res.update(pi_hom=float(a.mean()), pi_nonhom=float(b.mean()), n_hom=len(a), n_nonhom=len(b))
    # pairwise dS and upstream divergence
    ups = [up_seq.get(g, "") for g in ids]
    use_up = all(len(u) >= 40 for u in ups)
    aup = mafft(ups) if use_up else None
    dss, dns, ratios, dups = [], [], [], []
    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            a, b = nts[i], nts[j]
            keep = [k for k in range(0, len(a), 3) if a[k:k+3] != "---" and b[k:k+3] != "---"
                    and "*" not in str(Seq(a[k:k+3]).translate(table=11)) + str(Seq(b[k:k+3]).translate(table=11))]
            if len(keep) < 50:
                continue
            ca = CodonSeq("".join(a[k:k+3] for k in keep)); cb = CodonSeq("".join(b[k:k+3] for k in keep))
            try:
                dn, ds = cal_dn_ds(ca, cb, method="NG86")
            except Exception:
                continue
            if not (np.isfinite(ds) and ds >= 0):
                continue
            dss.append(ds); dns.append(dn)
            if use_up and 0.02 < ds < 1.5:
                du = jc(aup[i], aup[j])
                if np.isfinite(du):
                    dups.append(du); ratios.append(du / ds)
    res.update(dS_pair_median=float(np.median(dss)) if dss else np.nan,
               dN_pair_median=float(np.median(dns)) if dns else np.nan,
               up_len_median=float(np.median([len(u) for u in ups])),
               d_up_median=float(np.median(dups)) if len(dups) >= 3 else np.nan,
               up_constraint=float(np.median(ratios)) if len(ratios) >= 3 else np.nan)
    return res


if __name__ == "__main__":
    files = sorted(ALN.glob("*.fna"))
    with Pool(int(os.environ.get("THREADS", 2))) as p:
        rows = list(p.imap_unordered(analyse, files, chunksize=4))
    pd.DataFrame(rows).to_csv(OUT / "site_stats.tsv", sep="\t", index=False)
    print(len(rows), "OG-clade alignments analysed")
