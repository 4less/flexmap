#!/bin/bash
# Datasets for the flexmap benchmark: uppercase FASTA, one line per record, headers g<genome>_c<contig>.
set -euo pipefail
out=$1
norm() { awk -v g="$1" 'BEGIN{c=0} /^>/{if(c) printf "\n"; printf ">g%s_c%d\n", g, c++; next} {printf "%s", toupper($0)} END{printf "\n"}'; }
# D1: 450 random genomes of the GTDB downloader's raw set (mixed species)
ls ~/.gtdb_downloader/raw/*.fna.gz | shuf --random-source=<(yes) -n 450 > $out/d1_list.txt
i=0; while read f; do zcat "$f" | norm $i; i=$((i+1)); done < $out/d1_list.txt > $out/d1_diverse.fa
# D2: 300 Bacteroides uniformis genomes (one species: many cores shared by many genomes)
d=~/.gtdb_downloader/r232/genomes/taxonomy/Bacteria/Bacteroidota/Bacteroidia/Bacteroidales/Bacteroidaceae/Bacteroides/Bacteroides_uniformis
ls $d/*.fna* | shuf --random-source=<(yes) -n 300 > $out/d2_list.txt
i=0; while read f; do (zcat -f "$f") | norm $i; i=$((i+1)); done < $out/d2_list.txt > $out/d2_uniformis.fa
# D3: human GRCh37 chromosomes 1-8
awk '/^>/{p = ($1 ~ /^>NC_00000[1-8]\./)} p' ~/Databases/genomes/GRCh37/GCF_000001405.13_GRCh37_genomic.fna | norm 0 > $out/d3_human_chr1-8.fa
for f in $out/*.fa; do echo "$f $(grep -c '>' $f) records $(grep -v '>' $f | tr -d '\n' | wc -c) bases"; done
