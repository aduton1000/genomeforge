#!/usr/bin/env bash
# run_genomeforge.sh — launch (or resume) GenomeForge on a SLURM cluster. The Nextflow head
# job runs under sbatch so it survives logout; work/ and tmp live under the run directory
# (disk), never on a RAM-backed /tmp.
#
#   bash run_genomeforge.sh --input samplesheet.csv [extra nextflow args...]
#   bash run_genomeforge.sh --input samplesheet.csv -stub-run        # DAG dry run on the cluster
#
# Environment:
#   GENOMEFORGE_RUNDIR     run directory (default: $PWD/genomeforge_run)
#   GENOMEFORGE_PROFILE    nextflow profile (default: noguchi)
#   GENOMEFORGE_PARAMS     -params-file (default: none; the site profile carries the paths)
#   GENOMEFORGE_PARTITION  SLURM partition for the head job (default: global)
#
# Monitor:  tail -f $RUNDIR/nextflow.log      squeue -u $USER      $RUNDIR/results/pipeline_info/report.html
set -euo pipefail
PIPE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUNDIR="${GENOMEFORGE_RUNDIR:-$PWD/genomeforge_run}"
PROFILE="${GENOMEFORGE_PROFILE:-noguchi}"
PARAMS="${GENOMEFORGE_PARAMS:-}"
PARTITION="${GENOMEFORGE_PARTITION:-global}"
mkdir -p "$RUNDIR/work" "$RUNDIR/tmp" "$RUNDIR/logs"
[ -w "$RUNDIR" ] || { echo "run dir not writable: $RUNDIR" >&2; exit 1; }

STAMP=$(date +%Y%m%d_%H%M%S)
HEAD="$RUNDIR/logs/head_${STAMP}.sbatch"
cat > "$HEAD" <<EOF
#!/bin/bash
#SBATCH --job-name=genomeforge-head
#SBATCH --partition=${PARTITION}
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=${RUNDIR}/logs/head_${STAMP}.out
set -u
cd "${RUNDIR}"
export NXF_WORK="${RUNDIR}/work"
export TMPDIR="${RUNDIR}/tmp"
export NXF_TEMP="${RUNDIR}/tmp"
export NXF_OPTS="-Xms1g -Xmx6g -Djava.io.tmpdir=${RUNDIR}/tmp"
export NXF_ANSI_LOG=false
export NXF_SINGULARITY_CACHEDIR="\${NXF_SINGULARITY_CACHEDIR:-\$HOME/.genomeforge/container_cache}"
ulimit -n 65536 2>/dev/null || true
echo "GenomeForge head job \$SLURM_JOB_ID on \$(hostname) at \$(date)"
nextflow run "${PIPE}/main.nf" -profile ${PROFILE} ${PARAMS:+-params-file "${PARAMS}"} \\
    -resume -with-trace -work-dir "${RUNDIR}/work" $(printf '%q ' "$@") \\
    2>&1 | tee -a "${RUNDIR}/nextflow.log"
echo "head job finished at \$(date) with status \${PIPESTATUS[0]}"
EOF
echo "run dir : $RUNDIR"
echo "profile : $PROFILE ${PARAMS:+(params: $PARAMS)}"
echo "head job: $HEAD"
JOB=$(sbatch --parsable "$HEAD")
echo "submitted head job $JOB — follow with: tail -f $RUNDIR/nextflow.log"
