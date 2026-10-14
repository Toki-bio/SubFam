#!/usr/bin/env bash
# MeShClust2 sweeps on the three real sets, with the scorer of meshclust_sweep.py.
#   run_meshclust.sh MESHCLUST2_BIN OUTDIR B2.fa B2_CHUNKS.tsv B2_TRUTH.tsv L1.fa L1_CHUNKS.tsv L1_TRUTH.tsv ALU.fa ALU_CHUNKS.tsv ALU_TRUTH.tsv [THREADS]
# Inputs as in ../baseline_vsearch/run_sweeps.sh and ../sine_b2/make_b2.sh (chunks.tsv = SubFam's output on the same copies).
# MeShClust2: github.com/BioinformaticsToolsmith/MeShClust2 commit d7fdb42, built with g++ 13.3 after adding
# "-include cstdint -include limits -include algorithm" to CMAKE_CXX_FLAGS in its CMakeLists.txt.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); M=$1; O=$2; shift 2; T=${12:-8}; mkdir -p "$O"
python3 "$HERE/meshclust_sweep.py" "$M" "$1" "$2" "$3" rmsk --ids 0.60,0.70,0.80,0.90,0.95 --threads "$T" > "$O/b2_mesh.txt" 2>&1
python3 "$HERE/meshclust_sweep.py" "$M" "$4" "$5" "$6" rmsk --ids 0.80,0.90,0.95,0.98 --threads "$T" > "$O/l1_mesh.txt" 2>&1
python3 "$HERE/meshclust_sweep.py" "$M" "$7" "$8" "$9" Price --ids 0.90,0.95,0.98,0.99 --threads "$T" > "$O/alu_mesh.txt" 2>&1
