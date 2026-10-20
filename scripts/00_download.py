"""Download the input data from NCBI E-utilities (run once; needs internet).

  * chlorovirus genomes listed in data/accessions.txt  -> data/raw/chlorovirus_genomes.fna
  * metadata (length, organism, title)                  -> data/raw/chlorovirus_meta.tsv
  * green-algal proteomes (host + relatives)            -> data/raw/host_<species>.faa
  * PBCV-1 RefSeq protein annotation (gene names)       -> data/raw/PBCV1_annot.faa

Snapshot used in the project: NCBI, 2026-10-20 (query below).
Accession list was obtained with:  Chlorovirus[Organism] AND 250000:450000[SLEN]
(RefSeq NC_ duplicates of GenBank records removed).
"""
import json, time, urllib.parse, urllib.request
from pathlib import Path

E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
RAW = Path("data/raw"); RAW.mkdir(parents=True, exist_ok=True)
ALGAE = {
    "Chlorella_variabilis": '"Chlorella variabilis"[Organism]',
    "Micractinium_conductrix": '"Micractinium conductrix"[Organism]',
    "Coccomyxa_subellipsoidea": '"Coccomyxa subellipsoidea C-169"[Organism]',
    "Chlamydomonas_reinhardtii": '"Chlamydomonas reinhardtii"[Organism] AND refseq[filter]',
}


def get(url, data=None, tries=4):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, data=urllib.parse.urlencode(data).encode() if data else None)
            return urllib.request.urlopen(req, timeout=300).read().decode()
        except Exception:
            time.sleep(3 * (i + 1))
    raise RuntimeError(url)


accs = Path("data/accessions.txt").read_text().split()
with open(RAW / "chlorovirus_genomes.fna", "w") as fh:
    for i in range(0, len(accs), 20):
        fh.write(get(E + "efetch.fcgi", dict(db="nuccore", rettype="fasta", retmode="text", id=",".join(accs[i:i + 20]))))
        time.sleep(0.4)
summ = json.loads(get(E + "esummary.fcgi", dict(db="nuccore", retmode="json", id=",".join(accs))))["result"]
with open(RAW / "chlorovirus_meta.tsv", "w") as fh:
    fh.write("acc\tlen\torganism\tcreatedate\ttitle\n")
    for uid in summ["uids"]:
        r = summ[uid]
        fh.write(f"{r['accessionversion']}\t{r['slen']}\t{r['organism']}\t{r['createdate']}\t{r['title']}\n")
for sp, q in ALGAE.items():
    s = json.loads(get(E + "esearch.fcgi?" + urllib.parse.urlencode(dict(db="protein", retmode="json", usehistory="y", retmax=0, term=q))))["esearchresult"]
    with open(RAW / f"host_{sp}.faa", "w") as fh:
        for st in range(0, int(s["count"]), 5000):
            fh.write(get(E + "efetch.fcgi?" + urllib.parse.urlencode(dict(db="protein", rettype="fasta", retmode="text", query_key=s["querykey"], WebEnv=s["webenv"], retstart=st, retmax=5000))))
            time.sleep(0.4)
(RAW / "PBCV1_annot.faa").write_text(get(E + "efetch.fcgi?db=nuccore&rettype=fasta_cds_aa&retmode=text&id=NC_000852.5"))
print("done")
