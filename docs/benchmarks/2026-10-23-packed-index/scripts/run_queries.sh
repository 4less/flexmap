#!/bin/bash
# One process per index variant: load, then time lookups. Usage: run_queries.sh <dataset> <out.tsv>
set -euo pipefail
d=$1; out=$2
B=~/git/4less/flexmap-packed-index/target/release/packbench
I=/dev/shm/fmbench; Q=$(dirname $0)/../idx/$d.q
export FLEXMAP_QUERIES=${QUERIES:-15000000}; R=${ROUNDS:-2}
export FLEXMAP_LONG=$(dirname $0)/../idx/${d}_long.q
run() { echo "# $*" >&2; env "$@" | sed -u "s/^/$d\t/" >> $out; }
: > $out
run FLEXMAP_THP=0 $B query old heap $I/$d.old.blob $Q 1,8 $R        # flexalign's blob as it is (cebfd56)
run FLEXMAP_THP=1 $B query old heap $I/$d.old.blob $Q 1,8 $R        # + huge pages
run FLEXMAP_THP=1 $B query old mmap $I/$d.old.blob $Q 1,8 $R
run FLEXMAP_THP=0 $B query compact heap $I/$d.compact.pk $Q 1,8 $R
run FLEXMAP_THP=1 $B query compact heap $I/$d.compact.pk $Q 1,8 $R
run FLEXMAP_THP=1 $B query compact mmap $I/$d.compact.pk $Q 1,8 $R
if [ -f $I/$d.direct.pk ]; then run FLEXMAP_THP=1 $B query direct heap $I/$d.direct.pk $Q 1,8 $R; fi
