#!/usr/bin/env bash
# Full pipeline. Usage:  conda env create -f environment.yml && conda activate ncldv-gene-package && bash run_all.sh
# Runtime: ~1.5-2 h on 2 CPU cores (dominated by codeml in step 4).
set -euo pipefail
export THREADS=${THREADS:-2}
[ -s data/raw/chlorovirus_genomes.fna ] || python scripts/00_download.py   # NCBI snapshot (needs internet)
md5sum -c data/raw_checksums.md5 || echo "WARNING: raw data differ from the analysed snapshot"
mkdir -p work && python scripts/01_genes.py 2> work/genes.log
mkdir -p work/og
mmseqs easy-cluster work/genes/proteins.faa work/og/all work/og/tmp --min-seq-id 0.3 -c 0.7 --cov-mode 0 -s 7.5 --threads "$THREADS" > work/og/mmseqs.log
python scripts/02_species_tree.py
iqtree -s work/tree/concat.faa -m LG+F+G4 -fast -alrt 1000 -T "$THREADS" --prefix work/tree/species -seed 42 -redo > work/tree/iqtree.stdout
mkdir -p work/host
python scripts/utils_tag_algae.py > work/host/algae.faa
diamond makedb --in work/host/algae.faa -d work/host/algae --quiet
diamond blastp -q work/genes/proteins.faa -d work/host/algae -o work/host/vs_algae.tsv --sensitive -e 1e-5 -k 5 \
  --threads "$THREADS" --outfmt 6 qseqid sseqid pident length qstart qend qlen sstart send slen evalue bitscore --quiet
python scripts/03_classify.py
python scripts/04_dnds.py
python scripts/05_sites_upstream.py
python scripts/06_stats_figures.py
echo "Done: see results/figures, results/tables, results/stats.txt"
