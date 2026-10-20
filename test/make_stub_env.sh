#!/usr/bin/env bash
# Creates a fake reference/annotation/FASTQ layout so the pipeline DAG can be exercised with
# `-stub-run` anywhere (no real data, no containers). Prints the nextflow command to run.
set -euo pipefail
PIPE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
T="${1:-$(mktemp -d)}"
mkdir -p "$T/ref" "$T/db" "$T/vep/homo_sapiens/112_GRCh38" "$T/delivery/part1/01.RawData" "$T/delivery/part2/X-F004/01.RawData"
cd "$T"
touch ref/ref.fa ref/ref.fa.fai ref/ref.dict ref/ref.fa.0123 ref/ref.fa.amb ref/ref.fa.ann ref/ref.fa.bwt.2bit.64 ref/ref.fa.pac
printf 'chr1\t1000\nchr2\t1000\nchrX\t1000\nchrY\t500\nchrM\t100\n' > ref/ref.fa.fai
for v in dbsnp.chr.vcf.gz mills.vcf.gz known_indels.vcf.gz hapmap.vcf.gz omni.vcf.gz g1000.vcf.gz clinvar.vcf.gz sites.hg38.vcf.gz gnomad.exomes.v4.1.sites.chr1.vcf.bgz gnomad.exomes.v4.1.sites.chr22.vcf.bgz; do
  touch "db/$v" "db/$v.tbi"
done
mk() { d="$1/01.RawData/$2"; mkdir -p "$d"; lib="FKDN2603$((RANDOM % 900 + 100))-1A"; for u in "${@:3}"; do fc=${u%%:*}; ln=${u##*:}; for r in 1 2; do head -c 1200 /dev/urandom | gzip > "$d/$2_${lib}_${fc}_L${ln}_${r}.fq.gz"; done; done; }
mk delivery/part1 S001 23N5Y5LT4:5 23NF2CLT4:4
mk delivery/part1 S002 23N5Y3LT4:6
mk delivery/part1 S003 23N5Y5LT4:3 23NF2CLT4:3 23N2CNLT4:1
mk delivery/part2/X-F004 S004 252LYFLT4:8 253WYMLT4:7
# other conventions: bcl2fastq (lane in name, flowcell from header) and generic (both from header)
mkdir -p delivery/other
hdr() { printf '@A00551:123:%s:%s:1101:1000:2000 %s:N:0:ACGT\nACGT\n+\nIIII\n' "$1" "$2" "$3"; }
for r in 1 2; do hdr HXXXXDSXY 2 $r | gzip > "delivery/other/S005_S7_L002_R${r}_001.fastq.gz"; hdr HYYYYDSXY 3 $r | gzip > "delivery/other/S006.R${r}.fastq.gz"; done
printf 'sample_id,sex,phenotype\nS001,F,case\nS002,F,control\nS003,M,case\nS004,F,control\nS005,M,case\nS006,F,control\n' > meta.csv
python3 "$PIPE/bin/make_samplesheet.py" --roots delivery/part1 delivery/part2 delivery/other --metadata meta.csv --out samplesheet.csv --min-bytes 1
echo
echo "cd $T && nextflow run $PIPE/main.nf -profile stub -stub-run --input samplesheet.csv --outdir results --callforge_sif x --cohortqc_sif x --burden_sif x \\"
echo "  --fasta $T/ref/ref.fa --dbsnp $T/db/dbsnp.chr.vcf.gz --mills $T/db/mills.vcf.gz --known_indels $T/db/known_indels.vcf.gz \\"
echo "  --hapmap $T/db/hapmap.vcf.gz --omni $T/db/omni.vcf.gz --g1000 $T/db/g1000.vcf.gz --clinvar $T/db/clinvar.vcf.gz \\"
echo "  --gnomad_dir $T/db --somalier_sites $T/db/sites.hg38.vcf.gz --vep_cache $T/vep --chromosomes chr1,chr22,chrM"
