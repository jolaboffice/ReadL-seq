#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_DIR}"

POD5="data/example_lambda/example_lambda.pod5"
REF="data/example_lambda/lambda_forward.fasta"
LEVELS="data/resources/levels.txt"
OUTDIR="data/example_lambda/results"

mkdir -p "${OUTDIR}"

echo "========================================"
echo "ReadL-seq example pipeline"
echo "========================================"

echo "[1/3] Basecalling with Dorado..."
dorado basecaller \
  --disable-read-splitting \
  --no-trim \
  --emit-moves \
  dna_r10.4.1_e8.2_400bps_sup@v5.2.0 \
  "${POD5}" \
  > "${OUTDIR}/example.basecalled.bam"

echo "[2/3] Alignment with minimap2..."
samtools fastq -T "qs,du,ns,ts,mx,ch,st,TS,rn,fn,sm,sd,sv,dx,RG,mv" "${OUTDIR}/example.basecalled.bam" \
  | minimap2 -Y --secondary=no --end-bonus=0 -B 6 -O 8,48 -I 6G -E 2,1 -y -ax map-ont --MD "${REF}" - \
  | samtools view -b -F 0x900 \
  | samtools sort -o "${OUTDIR}/example.aligned.sorted.bam"

samtools index "${OUTDIR}/example.aligned.sorted.bam"

echo "[3/3] Running biotin_ssb_pipeline.py..."
python scripts/biotin_ssb_pipeline.py \
  --bam "${OUTDIR}/example.aligned.sorted.bam" \
  --pod5 "${POD5}" \
  --reference "${REF}" \
  --level-table "${LEVELS}" \
  --output "results"

echo "========================================"
echo "Done"
echo "Results saved in: ${OUTDIR}"
echo "========================================"
