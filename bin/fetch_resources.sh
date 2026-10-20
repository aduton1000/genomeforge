#!/usr/bin/env bash
# fetch_resources.sh — stage the hg38 resources the pipeline needs that are NOT yet in
# /hpc/opt/callforge/refs/annotation_db: VQSR training sets + GATK WGS calling regions.
# Idempotent, resumable (curl -C -), verifies gzip integrity. ~2.5 GB total.
#
#   bash bin/fetch_resources.sh [DEST_DIR]      (default: /hpc/opt/callforge/refs/annotation_db)
set -euo pipefail
DEST="${1:-/hpc/opt/callforge/refs/annotation_db}"
BASE="https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0"
FILES=(
  hapmap_3.3.hg38.vcf.gz hapmap_3.3.hg38.vcf.gz.tbi
  1000G_omni2.5.hg38.vcf.gz 1000G_omni2.5.hg38.vcf.gz.tbi
  1000G_phase1.snps.high_confidence.hg38.vcf.gz 1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
  wgs_calling_regions.hg38.interval_list
)
mkdir -p "$DEST"; cd "$DEST"
for f in "${FILES[@]}"; do
  if [ -s "$f" ] && { [[ "$f" != *.gz ]] || gzip -t "$f" 2>/dev/null; }; then echo "have $f"; continue; fi
  echo "fetching $f"
  curl -fL --retry 5 --retry-delay 10 -C - -o "$f.part" "$BASE/$f"
  if [[ "$f" == *.gz ]]; then gzip -t "$f.part"; fi
  mv "$f.part" "$f"
done
echo
echo "contig style check (must all be chr-style):"
for f in hapmap_3.3.hg38.vcf.gz 1000G_omni2.5.hg38.vcf.gz 1000G_phase1.snps.high_confidence.hg38.vcf.gz; do
  printf '  %-50s %s\n' "$f" "$(tabix -l "$f" 2>/dev/null | head -1 || zcat "$f" | grep -v '^#' | head -1 | cut -f1)"
done
printf '  %-50s %s intervals\n' wgs_calling_regions.hg38.interval_list "$(grep -vc '^@' wgs_calling_regions.hg38.interval_list)"
ls -la "$DEST"/{hapmap_3.3,1000G_omni2.5,1000G_phase1.snps}*.vcf.gz* "$DEST"/wgs_calling_regions.hg38.interval_list
echo "done -> $DEST"
