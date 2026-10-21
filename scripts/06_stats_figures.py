"""Statistics and publication figures.

Inputs : results/tables/{genome_clades,orthogroups_classified}.tsv, work/dnds/dnds_raw.tsv,
         work/sites/site_stats.tsv, work/genes_annotated.tsv, work/tree/species.treefile
Outputs: results/figures/Fig{1..3}.{pdf,png}, results/tables/*.tsv, results/stats.txt
"""
import re
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from Bio import Phylo

FIG = Path("results/figures"); FIG.mkdir(parents=True, exist_ok=True)
TAB = Path("results/tables")
LOG = open("results/stats.txt", "w")


def say(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n")


# ---------- style ----------
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 9,
                     "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                     "legend.frameon": False, "pdf.fonttype": 42, "svg.fonttype": "none"})
CLS = ["Core viral", "Core algal-like", "Accessory algal-like", "Accessory ORFan"]
CLS_LAB = ["Core\nviral", "Core\nalgal", "Acc.\nalgal\n(stolen)", "Acc.\nORFan"]
COL = {"Core viral": "#2a78d6", "Accessory algal-like": "#eb6834", "Core algal-like": "#1baf7a", "Accessory ORFan": "#eda100"}
CLADES = ["NC64A", "Osy", "SAG", "Pbi"]
CCOL = {"NC64A": "#e87ba4", "Osy": "#008300", "SAG": "#6250d6", "Pbi": "#e34948", "unassigned": "#8a8986"}
MARK = {"NC64A": "o", "Osy": "s", "SAG": "^", "Pbi": "D"}
INK2 = "#52514e"


def save(fig, name):
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def cliffs_delta(a, b):
    a, b = np.asarray(a), np.asarray(b)
    gt = (a[:, None] > b[None, :]).sum(); lt = (a[:, None] < b[None, :]).sum()
    return (gt - lt) / (len(a) * len(b))


def pstr(p):
    return "p < 1e-10" if p < 1e-10 else (f"p = {p:.1e}" if p < 0.001 else f"p = {p:.3f}")


# ---------- data ----------
cl = pd.read_csv(TAB / "genome_clades.tsv", sep="\t")
ogt = pd.read_csv(TAB / "orthogroups_classified.tsv", sep="\t")
dn = pd.read_csv("work/dnds/dnds_raw.tsv", sep="\t")
st = pd.read_csv("work/sites/site_stats.tsv", sep="\t")
ga = pd.read_csv("work/genes_annotated.tsv", sep="\t")

d = dn[dn.error.isna()].merge(st.drop(columns=["n_seq"]), on=["clade", "og"], how="left").merge(
    ogt[["og", "class", "prevalence", "algal_frac", "Cv_pident_median", "algal_pident_median", "aa_len_median"]], on="og")
n0 = len(d)
d = d[(d.dS_pair_median > 0.01) & (d.dS_pair_median < 1.5) & (d.omega < 10)].copy()
say(f"# dN/dS dataset: {n0} clade-OG estimates, {len(d)} after dS filter (0.01 < median pairwise dS < 1.5)")
d["log_omega"] = np.log10(d.omega.clip(lower=1e-3))
d["LRT"] = 2 * (d.lnL - d.lnL_w1)
d["p_LRT"] = stats.chi2.sf(d.LRT.clip(lower=0), 1)
d["q_LRT"] = multipletests(d.p_LRT, method="fdr_bh")[1]
d["regime"] = np.where(d.q_LRT >= 0.05, "not different from neutral",
                       np.where(d.omega < 1, "purifying (omega < 1)", "omega > 1"))
d.to_csv(TAB / "dnds_per_gene.tsv", sep="\t", index=False)

# ---------- descriptive ----------
say("\n## Genomes per clade\n" + cl.clade.value_counts().to_string())
say("\n## Orthogroups per class (global)\n" + ogt["class"].value_counts().to_string())
gpc = ga.merge(ogt[["og", "class"]], on="og").groupby(["genome", "class"]).size().unstack(fill_value=0)
say("\n## Genes per genome by class (median [min-max])")
for c in CLS:
    say(f"  {c}: {gpc[c].median():.0f} [{gpc[c].min()}-{gpc[c].max()}]  ({(gpc[c] / gpc.sum(1)).median() * 100:.1f}% of genes)")
say("\n## Estimates per clade x class\n" + d.groupby(["clade", "class"]).size().unstack().to_string())

# ---------- H1: are stolen genes neutral? ----------
say("\n# H1. Do accessory algal-like ('stolen') genes evolve neutrally (omega = 1)?")
summ = d.groupby("class").agg(n=("omega", "size"), median_omega=("omega", "median"),
                              q25=("omega", lambda x: x.quantile(.25)), q75=("omega", lambda x: x.quantile(.75)))
say(summ.loc[CLS].round(3).to_string())
reg = pd.crosstab(d["class"], d.regime).loc[CLS]
say("\nPer-gene LRT (codeml M0 vs omega=1, BH q<0.05):\n" + reg.to_string())
say((reg.div(reg.sum(1), axis=0) * 100).round(1).to_string())
for c in CLS:
    x = d.loc[d["class"] == c, "log_omega"]
    w = stats.wilcoxon(x)  # H0: median log10(omega) = 0, i.e. omega = 1
    say(f"  one-sample Wilcoxon log10(omega)=0, {c}: n={len(x)}, {pstr(w.pvalue)}")

say("\n## Class differences within each clade (Kruskal-Wallis; MWU stolen vs others, Cliff's delta)")
pw = []
for cd in CLADES:
    sub = d[d.clade == cd]
    kw = stats.kruskal(*[sub.loc[sub["class"] == c, "omega"] for c in CLS])
    say(f"  {cd}: KW H={kw.statistic:.1f}, {pstr(kw.pvalue)}; medians " +
        ", ".join(f"{c}={sub.loc[sub['class'] == c, 'omega'].median():.3f}" for c in CLS))
    a = sub.loc[sub["class"] == "Accessory algal-like", "omega"]
    for c in ["Core viral", "Core algal-like", "Accessory ORFan"]:
        b = sub.loc[sub["class"] == c, "omega"]
        u = stats.mannwhitneyu(a, b)
        pw.append(dict(clade=cd, comparison=f"Accessory algal-like vs {c}", n1=len(a), n2=len(b),
                       median1=a.median(), median2=b.median(), U=u.statistic, p=u.pvalue, cliffs_delta=cliffs_delta(a, b)))
pw = pd.DataFrame(pw); pw["q"] = multipletests(pw.p, method="fdr_bh")[1]
pw.to_csv(TAB / "omega_pairwise_tests.tsv", sep="\t", index=False)
say(pw.round(4).to_string(index=False))

say("\n## Pooled model: log10(omega) ~ class + clade + log10(dS) + log10(length); SE clustered by OG")
d["log_dS"] = np.log10(d.dS_pair_median); d["log_len"] = np.log10(d.n_codons)
d["gclass"] = d["class"]
m = smf.ols("log_omega ~ C(gclass, Treatment('Core viral')) + C(clade) + log_dS + log_len", data=d).fit(
    cov_type="cluster", cov_kwds={"groups": pd.factorize(d.og)[0]})
say(m.summary().tables[1].as_text())
say(f"  R2 = {m.rsquared:.3f}, n = {int(m.nobs)}")

say("\n## Recency of acquisition: omega vs identity to the Chlorella variabilis homolog (algal-like genes)")
alg = d[d["class"].str.contains("algal") & d.Cv_pident_median.notna()]
for c in ["Accessory algal-like", "Core algal-like"]:
    s = alg[alg["class"] == c]
    r = stats.spearmanr(s.Cv_pident_median, s.omega)
    say(f"  {c}: Spearman rho = {r.statistic:.3f}, {pstr(r.pvalue)}, n = {len(s)}")
r_all = stats.spearmanr(alg.Cv_pident_median, alg.omega)
say(f"  all algal-like: rho = {r_all.statistic:.3f}, {pstr(r_all.pvalue)}, n = {len(alg)}")

# ---------- H2: conserved regions ----------
say("\n# H2. Are there conserved regions inside stolen genes?")
h = d[d.p_hetero.notna()].copy()
h["q_hetero"] = multipletests(h.p_hetero, method="fdr_bh")[1]
het = h.groupby("class").agg(n=("q_hetero", "size"), frac_heterogeneous=("q_hetero", lambda x: (x < 0.05).mean()),
                             median_ratio=("hetero_ratio", "median"))
say("Conservation heterogeneity along the protein (window variance vs permutation; BH q<0.05):\n" + het.loc[CLS].round(3).to_string())
ct = pd.crosstab(h["class"], h.q_hetero < 0.05).loc[CLS]
chi = stats.chi2_contingency(ct)
say(f"  chi2 = {chi[0]:.1f}, df = {chi[2]}, {pstr(chi[1])}")
hh = d[d.pi_hom.notna() & d.pi_nonhom.notna()].copy()
hh["ratio"] = (hh.pi_nonhom + 1e-3) / (hh.pi_hom + 1e-3)
for c in ["Accessory algal-like", "Core algal-like"]:
    s = hh[hh["class"] == c]
    if len(s) >= 5:
        w = stats.wilcoxon(s.pi_hom, s.pi_nonhom)
        say(f"  {c}: pi_aa host-homologous = {s.pi_hom.median():.4f} vs non-homologous = {s.pi_nonhom.median():.4f}; "
            f"paired Wilcoxon {pstr(w.pvalue)}, n = {len(s)}; non-hom more variable in {(s.pi_nonhom > s.pi_hom).mean() * 100:.0f}%")
hh.to_csv(TAB / "homologous_region_conservation.tsv", sep="\t", index=False)

say("\n# H3. Is the upstream (promoter) region of stolen genes constrained like that of viral genes?")
u = d[d.up_constraint.notna()]
for c in CLS:
    x = u.loc[u["class"] == c, "up_constraint"]
    if len(x) >= 5:
        w = stats.wilcoxon(np.log(x))
        say(f"  {c}: n = {len(x)}, median d_up/dS = {x.median():.3f} [IQR {x.quantile(.25):.2f}-{x.quantile(.75):.2f}], "
            f"Wilcoxon vs 1: {pstr(w.pvalue)}; median upstream length {u.loc[u['class'] == c, 'up_len_median'].median():.0f} nt")
kw = stats.kruskal(*[u.loc[u["class"] == c, "up_constraint"] for c in CLS])
say(f"  KW across classes: H = {kw.statistic:.1f}, {pstr(kw.pvalue)}")
a = u.loc[u["class"] == "Accessory algal-like", "up_constraint"]; b = u.loc[u["class"] == "Core viral", "up_constraint"]
mw = stats.mannwhitneyu(a, b)
say(f"  stolen vs core viral: MWU {pstr(mw.pvalue)}, Cliff's delta = {cliffs_delta(a, b):.2f}")

# ---------- decay: truncated copies ----------
say("\n# Gene decay: fraction of truncated copies (<80% of OG median length) per clade-OG")
gg = ga.merge(ogt[["og", "class"]], on="og")
med = gg.groupby(["clade", "og"]).aa_len.transform("median")
gg["trunc"] = gg.aa_len < 0.8 * med
tr = gg[gg.clade.isin(CLADES)].groupby(["clade", "og", "class"]).agg(n=("trunc", "size"), trunc=("trunc", "mean")).reset_index()
tr = tr[tr.n >= 6]
say(tr.groupby("class").trunc.agg(["size", "mean", "median"]).loc[CLS].round(3).to_string())
kw = stats.kruskal(*[tr.loc[tr["class"] == c, "trunc"] for c in CLS])
say(f"  KW H = {kw.statistic:.1f}, {pstr(kw.pvalue)}")

# ---------- named examples: PBCV-1 (JF411744.1 = RefSeq NC_000852.5) annotations ----------
ann = []
for line in open("data/raw/PBCV1_annot.faa"):
    if line.startswith(">"):
        loc = re.search(r"\[location=([^\]]+)\]", line).group(1)
        prot = re.search(r"\[protein=([^\]]+)\]", line); tag = re.search(r"\[locus_tag=([^\]]+)\]", line)
        nums = [int(x) for x in re.findall(r"\d+", loc)]
        strand = "-" if loc.startswith("complement") else "+"
        stop = min(nums) if strand == "-" else max(nums)
        ann.append(dict(strand=strand, stop=stop, product=prot.group(1) if prot else "", locus_tag=tag.group(1) if tag else ""))
ann = pd.DataFrame(ann)
pb = ga[ga.genome == "JF411744.1"].copy()
pb["stop"] = np.where(pb.strand == "-", pb.start, pb.end)
pb = pb.merge(ann, on=["strand", "stop"], how="left")
named = pb[["og", "gene", "locus_tag", "product"]].merge(
    d[d.clade == "NC64A"][["og", "class", "omega", "q_LRT", "Cv_pident_median", "pi_hom", "pi_nonhom", "up_constraint"]], on="og")
named = named[named["class"] == "Accessory algal-like"].drop_duplicates("gene").sort_values("omega")
named.to_csv(TAB / "stolen_genes_PBCV1_named.tsv", sep="\t", index=False)
with open("results/stats.txt", "a") as fh:
    fh.write("\n# Accessory algal-like genes present in PBCV-1 (NC64A clade estimates)\n")
    fh.write(named.round(3).to_string(index=False) + "\n")
print(named.round(3).to_string(index=False))

# ===================== FIGURE 1: data & phylogeny =====================
tree = Phylo.read("work/tree/species.treefile", "newick"); tree.root_at_midpoint(); tree.ladderize()
clade_of = dict(zip(cl.genome, cl.clade))
fig = plt.figure(figsize=(7.2, 4.2))
ax = fig.add_axes([0.0, 0.02, 0.42, 0.93])
depths = tree.depths(); term = tree.get_terminals()
ypos = {t: i for i, t in enumerate(term)}
def setY(c):
    if c.is_terminal(): return ypos[c]
    ypos[c] = np.mean([setY(ch) for ch in c.clades]); return ypos[c]
setY(tree.root)
for c in tree.find_clades():
    x = depths[c]
    for ch in c.clades:
        ax.plot([x, x], [ypos[c], ypos[ch]], color="#3a3a38", lw=0.5)
        ax.plot([x, depths[ch]], [ypos[ch], ypos[ch]], color="#3a3a38", lw=0.5)
    if not c.is_terminal() and c.confidence is not None and c.confidence >= 80 and len(c.get_terminals()) > 8:
        ax.plot(x, ypos[c], "o", ms=2.6, color="#3a3a38")
xmax = max(depths[t] for t in term)
for t in term:
    ax.plot(depths[t], ypos[t], "o", ms=2.4, color=CCOL[clade_of.get(t.name, "unassigned")], mec="none")
for cd in CLADES:
    ys = [ypos[t] for t in term if clade_of.get(t.name) == cd]
    ax.plot([xmax * 1.04] * 2, [min(ys), max(ys)], color=CCOL[cd], lw=3, solid_capstyle="butt")
    host = {"NC64A": "C. variabilis NC64A", "Osy": "C. variabilis Syngen", "SAG": "C. heliozoae SAG", "Pbi": "M. conductrix Pbi"}[cd]
    ax.text(xmax * 1.07, np.mean(ys), f"{cd}-viruses (n={len(ys)})\nhost: {host}", va="center", fontsize=6.5, color="#0b0b0b")
ax.plot([0, 0.1], [-2, -2], color="#3a3a38", lw=0.8); ax.text(0.05, -3, "0.1 subst./site", ha="center", va="bottom", fontsize=6)
ax.set_ylim(len(term) + 1, -7); ax.set_xlim(-0.02, xmax * 1.75); ax.axis("off")
ax.set_title("A   Chlorovirus phylogenomic tree", loc="left", fontweight="bold")
ax.text(0, len(term) + 0.5, f"{len(term)} genomes, ML on {pd.read_csv('work/tree/markers.txt', header=None).shape[0]} single-copy proteins\n(LG+F+G4); dots: SH-aLRT >= 80", fontsize=6, color=INK2, va="top")

ax = fig.add_axes([0.53, 0.58, 0.45, 0.33])
bins = np.linspace(0, 1, 21)
bottom = np.zeros(20)
for c in CLS:
    hcnt, _ = np.histogram(ogt.loc[ogt["class"] == c, "prevalence"], bins=bins)
    ax.bar(bins[:-1], hcnt, width=0.05 * 0.92, bottom=bottom, align="edge", color=COL[c], label=c, edgecolor="white", lw=0.4)
    bottom += hcnt
YMAX = 200
ax.set_ylim(0, YMAX); ax.set_xlabel("Prevalence across 132 genomes"); ax.set_ylabel("Orthogroups")
if bottom[0] > YMAX:
    ax.text(0.025, YMAX * 0.97, f"{int(bottom[0])} ↑", ha="center", va="top", fontsize=6, color="#0b0b0b")
ax.axvline(0.9, color=INK2, lw=0.6, ls="--"); ax.text(0.885, YMAX * 0.95, "core\n≥ 0.9", ha="right", va="top", fontsize=6, color=INK2)
ax.set_title("B   Pan-genome: 1,642 orthogroups", loc="left", fontweight="bold")
ax.legend(fontsize=6, loc="upper center", ncol=2, bbox_to_anchor=(0.5, -0.32), handlelength=1)

ax = fig.add_axes([0.53, 0.08, 0.45, 0.24])
comp = gpc[CLS].div(gpc.sum(1), axis=0).join(cl.set_index("genome").clade)
comp = comp[comp.clade.isin(CLADES)].groupby("clade")[CLS].mean().loc[CLADES]
left = np.zeros(len(comp))
for c in CLS:
    ax.barh(range(len(comp)), comp[c], left=left, color=COL[c], edgecolor="white", lw=0.6, height=0.7)
    for i, v in enumerate(comp[c]):
        if v > 0.05:
            ax.text(left[i] + v / 2, i, f"{v * 100:.0f}%", ha="center", va="center", fontsize=6, color="white" if c == "Core viral" else "#0b0b0b")
    left += comp[c].values
ax.set_yticks(range(len(comp))); ax.set_yticklabels([f"{c}" for c in comp.index]); ax.invert_yaxis()
ax.set_xlim(0, 1); ax.set_xlabel("Mean fraction of genes per genome")
ax.set_title("C   Gene content by class", loc="left", fontweight="bold")
save(fig, "Fig1_data_phylogeny")

# ===================== FIGURE 2: selection =====================
fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.8), gridspec_kw=dict(width_ratios=[1.3, 1, 1], wspace=0.5))
ax = axs[0]
rng = np.random.default_rng(1)
for i, c in enumerate(CLS):
    x = d.loc[d["class"] == c, "omega"].values
    ax.scatter(i + rng.uniform(-0.28, 0.28, len(x)), x, s=3, color=COL[c], alpha=0.35, lw=0, rasterized=True)
    q1, q2, q3 = np.quantile(x, [.25, .5, .75])
    ax.plot([i - 0.32, i + 0.32], [q2, q2], color="#0b0b0b", lw=1.6)
    ax.plot([i, i], [q1, q3], color="#0b0b0b", lw=0.9)
    ax.text(i, 2.2, f"{q2:.2f}\nn={len(x)}", ha="center", va="bottom", fontsize=6)
ax.axhline(1, color="#c0392b", lw=0.8, ls="--"); ax.text(-0.45, 1.08, "neutral (ω = 1)", fontsize=6, color="#c0392b", ha="left", va="bottom")
ax.set_yscale("log"); ax.set_ylim(5e-3, 8); ax.set_xticks(range(4)); ax.set_xticklabels(CLS_LAB, fontsize=6.5)
ax.set_ylabel("dN/dS (ω, codeml M0)")
ax.set_title("A   dN/dS by gene class", loc="left", fontweight="bold")

ax = axs[1]
for cd in CLADES:
    sub = d[d.clade == cd]
    med = [sub.loc[sub["class"] == c, "omega"].median() for c in CLS]
    ax.plot(range(4), med, color=CCOL[cd], lw=1, marker=MARK[cd], ms=4, label=f"{cd}")
ax.axhline(1, color="#c0392b", lw=0.8, ls="--")
ax.set_yscale("log"); ax.set_ylim(0.015, 1.5); ax.set_xticks(range(4)); ax.set_xticklabels(CLS_LAB, fontsize=6.5)
ax.set_ylabel("median ω per clade")
ax.legend(fontsize=6, title="clade", title_fontsize=6, loc="lower left", ncol=2)
ax.set_title("B   Four clades", loc="left", fontweight="bold")

ax = axs[2]
RG = ["purifying (omega < 1)", "not different from neutral", "omega > 1"]
RC = {"purifying (omega < 1)": "#2a78d6", "not different from neutral": "#c3c2b7", "omega > 1": "#e34948"}
frac = pd.crosstab(d["class"], d.regime, normalize="index").reindex(columns=RG, fill_value=0).loc[CLS]
bot = np.zeros(4)
for r in RG:
    ax.bar(range(4), frac[r], bottom=bot, color=RC[r], width=0.7, edgecolor="white", lw=0.6,
           label={"purifying (omega < 1)": "ω < 1 (purifying)", "not different from neutral": "n.s. (≈ neutral)", "omega > 1": "ω > 1"}[r])
    for i, v in enumerate(frac[r]):
        if v > 0.06:
            ax.text(i, bot[i] + v / 2, f"{v * 100:.0f}%", ha="center", va="center", fontsize=6, color="white" if r != "not different from neutral" else "#0b0b0b")
    bot += frac[r].values
ax.set_xticks(range(4)); ax.set_xticklabels(CLS_LAB, fontsize=6.5); ax.set_ylim(0, 1); ax.set_ylabel("Fraction of genes")
ax.legend(fontsize=5.8, loc="upper center", bbox_to_anchor=(0.45, -0.27), ncol=1, handlelength=1)
ax.set_title("C   LRT vs ω = 1", loc="left", fontweight="bold")
save(fig, "Fig2_selection")

# ===================== FIGURE 3: regions =====================
fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.7), gridspec_kw=dict(wspace=0.45, width_ratios=[0.9, 1.3, 1.1]))
ax = axs[0]
s = hh[hh["class"].str.contains("algal")]
for _, r in s.iterrows():
    ax.plot([0, 1], [r.pi_hom + 1e-3, r.pi_nonhom + 1e-3], color=COL[r["class"]], lw=0.5, alpha=0.5)
for i, col in enumerate(["pi_hom", "pi_nonhom"]):
    ax.plot([i - 0.15, i + 0.15], [s[col].median() + 1e-3] * 2, color="#0b0b0b", lw=2)
ax.set_yscale("log"); ax.set_xticks([0, 1]); ax.set_xticklabels(["homologous\nto algal protein", "viral-specific\nregion"], fontsize=6.5)
ax.set_xlim(-0.4, 1.4); ax.set_ylabel("aa diversity per site (π$_{aa}$)")
pp = stats.wilcoxon(s.pi_hom, s.pi_nonhom).pvalue
ax.set_title("A   Host-derived core is\n     more conserved", loc="left", fontweight="bold")
ax.text(0.5, 0.97, f"paired Wilcoxon {pstr(pp)}\nn = {len(s)} gene-clade pairs", transform=ax.transAxes, ha="center", va="top", fontsize=6)
ax.legend(handles=[Line2D([], [], color=COL[c], lw=1.2, label=c) for c in ["Core algal-like", "Accessory algal-like"]],
          fontsize=5.8, loc="upper center", bbox_to_anchor=(0.5, -0.2))

ax = axs[1]
# example: the stolen gene with strongest heterogeneity in NC64A
ex = h[(h["class"] == "Accessory algal-like")].sort_values(["q_hetero", "hetero_ratio"], ascending=[True, False])
ex = ex[ex.clade == "NC64A"].iloc[0] if (ex.clade == "NC64A").any() else ex.iloc[0]
from Bio import SeqIO
from Bio.Seq import Seq
recs = list(SeqIO.parse(f"work/dnds/aln/{ex.clade}_{ex.og}.fna", "fasta"))
prots = ["".join("-" if str(r.seq)[i:i + 3] == "---" else str(Seq(str(r.seq)[i:i + 3]).translate(table=11)) for i in range(0, len(r.seq), 3)) for r in recs]
L = len(prots[0]); pis = []; cov = np.zeros(L)
sub_g = ga.set_index("gene").loc[[r.id for r in recs]]
for i in range(L):
    aa = [p[i] for p in prots if p[i] not in "-X*"]
    if len(aa) < 0.8 * len(prots): pis.append(np.nan); continue
    vc = pd.Series(aa).value_counts().values; m_ = len(aa)
    pis.append(1 - (vc * (vc - 1)).sum() / (m_ * (m_ - 1)))
for p, (gid, row) in zip(prots, sub_g.iterrows()):
    if not row.algal_hit: continue
    k = 0
    for i, a in enumerate(p):
        if a != "-":
            k += 1
            if row.qstart <= k <= row.qend: cov[i] += 1
pis = np.array(pis); xs = np.arange(L)
smooth = pd.Series(pis).rolling(15, center=True, min_periods=8).mean()
hom = cov >= 0.5 * sub_g.algal_hit.sum()
ax.fill_between(xs, 0, 1, where=hom, transform=ax.get_xaxis_transform(), color="#1baf7a", alpha=0.13, lw=0)
ax.bar(xs, pis, width=1, color="#c3c2b7")
ax.plot(xs, smooth, color="#eb6834", lw=1.4)
ax.set_xlim(0, L); ax.set_xlabel("Alignment position (codon)"); ax.set_ylabel("π$_{aa}$")
pb = ga[ga.og == ex.og]
lt = named.loc[named.og == ex.og, "locus_tag"]
nm = lt.iloc[0].replace("PBCV1_", "PBCV-1 ") if len(lt) and isinstance(lt.iloc[0], str) else ex.og
ax.set_title(f"B   Example: {nm}", loc="left", fontweight="bold")
ax.set_ylim(0, 1.0); ax.set_yticks([0, 0.2, 0.4, 0.6])
ax.text(0.98, 0.98, f"{ex.clade}: ω = {ex.omega:.2f}; {ex.Cv_pident_median:.0f}% id. to C. variabilis\nshaded: homologous to algal protein\nline: 15-codon mean (q = {ex.q_hetero:.0e})",
        transform=ax.transAxes, va="top", ha="right", fontsize=5.6, color=INK2, bbox=dict(fc="white", ec="none", alpha=0.9, pad=1))

ax = axs[2]
for i, c in enumerate(CLS):
    x = u.loc[u["class"] == c, "up_constraint"].values
    ax.scatter(i + rng.uniform(-0.28, 0.28, len(x)), x, s=3, color=COL[c], alpha=0.4, lw=0, rasterized=True)
    q1, q2, q3 = np.quantile(x, [.25, .5, .75])
    ax.plot([i - 0.32, i + 0.32], [q2, q2], color="#0b0b0b", lw=1.6); ax.plot([i, i], [q1, q3], color="#0b0b0b", lw=0.9)
    ax.text(i, 12, f"{q2:.2f}\nn={len(x)}", ha="center", va="bottom", fontsize=6)
ax.axhline(1, color="#c0392b", lw=0.8, ls="--")
ax.text(-0.45, 1.1, "= dS (neutral)", fontsize=5.8, color="#c0392b", ha="left", va="bottom")
ax.set_yscale("log"); ax.set_ylim(0.02, 40); ax.set_xticks(range(4)); ax.set_xticklabels(CLS_LAB, fontsize=6.5)
ax.set_ylabel("upstream divergence / dS")
ax.set_title("C   Upstream constraint", loc="left", fontweight="bold")
save(fig, "Fig3_regions_promoters")

LOG.close()
