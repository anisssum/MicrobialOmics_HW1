"""Define host-range clades on the species tree and classify orthogroups.

Clades: the largest monophyletic group containing exactly one reference virus
  NC64A-viruses  (PBCV-1, JF411744.1)       host: Chlorella variabilis NC64A
  SAG-viruses    (ATCV-1, EF101928.1)       host: Chlorella heliozoae SAG 3.83
  Pbi-viruses    (FR483,  DQ890022.1)       host: Micractinium conductrix Pbi
  Osy-viruses    (OSy-NE5, KX857749.1)      host: Chlorella variabilis Syngen 2-3

Orthogroup classes (per global OG):
  prevalence = fraction of all genomes carrying the OG
  algal-like = >=50% of members have a DIAMOND hit (e<=1e-5) in green-algal proteomes
  core     : prevalence >= 0.9      accessory : prevalence < 0.9
  -> Core viral | Core algal-like | Accessory algal-like ("stolen") | Accessory ORFan
"""
import collections
from pathlib import Path
import pandas as pd
from Bio import Phylo

REFS = {"NC64A": "JF411744.1", "SAG": "EF101928.1", "Pbi": "DQ890022.1", "Osy": "KX857749.1"}
out = Path("results/tables"); out.mkdir(parents=True, exist_ok=True)

tree = Phylo.read("work/tree/species.treefile", "newick")
tree.root_at_midpoint()
leaves = [t.name for t in tree.get_terminals()]
assign = {}
for clade_name, ref in REFS.items():
    best = None
    for node in tree.find_clades(order="level"):
        names = {t.name for t in node.get_terminals()}
        if ref in names and not any(r in names for k, r in REFS.items() if k != clade_name):
            if best is None or len(names) > len(best):
                best = names
    for n in best:
        assign[n] = clade_name
clades = pd.DataFrame({"genome": leaves, "clade": [assign.get(l, "unassigned") for l in leaves]})
meta = pd.read_csv("data/raw/chlorovirus_meta.tsv", sep="\t").rename(columns={"acc": "genome"})
clades = clades.merge(meta[["genome", "len", "organism", "title"]], on="genome", how="left")
clades.to_csv(out / "genome_clades.tsv", sep="\t", index=False)
print(clades.clade.value_counts())

genes = pd.read_csv("work/genes/genes.tsv", sep="\t")
og = pd.read_csv("work/og/orthogroups.tsv", sep="\t")
genes = genes.merge(og, on="gene").merge(clades[["genome", "clade"]], on="genome")

hits = pd.read_csv("work/host/vs_algae.tsv", sep="\t", header=None,
                   names="qseqid sseqid pident length qstart qend qlen sstart send slen evalue bitscore".split())
hits["species"] = hits.sseqid.str.split("|").str[0]
best = hits.sort_values("bitscore", ascending=False).drop_duplicates("qseqid")
bestCv = hits[hits.species == "Chlorella_variabilis"].sort_values("bitscore", ascending=False).drop_duplicates("qseqid")
genes["algal_hit"] = genes.gene.isin(best.qseqid)
genes = genes.merge(best[["qseqid", "pident", "qstart", "qend", "qlen", "species"]].rename(
    columns={"qseqid": "gene", "pident": "algal_pident", "species": "algal_best_sp"}), on="gene", how="left")
genes = genes.merge(bestCv[["qseqid", "pident"]].rename(columns={"qseqid": "gene", "pident": "Cv_pident"}), on="gene", how="left")
genes.to_csv("work/genes_annotated.tsv", sep="\t", index=False)

n_gen = clades.genome.nunique()
rows = []
for o, d in genes.groupby("og"):
    prev = d.genome.nunique() / n_gen
    alg = d.algal_hit.mean()
    rows.append(dict(og=o, n_genes=len(d), n_genomes=d.genome.nunique(), prevalence=prev,
                     algal_frac=alg, algal_pident_median=d.algal_pident.median(),
                     Cv_pident_median=d.Cv_pident.median(),
                     aa_len_median=d.aa_len.median(),
                     **{f"prev_{c}": d[d.clade == c].genome.nunique() / max(1, (clades.clade == c).sum()) for c in REFS}))
ogt = pd.DataFrame(rows)
core = ogt.prevalence >= 0.9
algal = ogt.algal_frac >= 0.5
ogt["class"] = "Accessory ORFan"
ogt.loc[core & ~algal, "class"] = "Core viral"
ogt.loc[core & algal, "class"] = "Core algal-like"
ogt.loc[~core & algal, "class"] = "Accessory algal-like"
ogt.to_csv(out / "orthogroups_classified.tsv", sep="\t", index=False)
print(ogt["class"].value_counts())
print(ogt.groupby("class").n_genes.sum())
