#!/bin/bash

# Get the current dir.
if [ -n "$BASH_VERSION" ]; then
    DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
elif [ -n "$ZSH_VERSION" ]; then
    DIR=${0:a:h}  # https://unix.stackexchange.com/a/115431
else
	echo "Error: Unknown shell; cannot determine path to integration-mapping-pipeline"
fi

export INTEGRATION_MAPPING_REPO_DIR="${DIR}"

# ── Run the Snakemake pipeline ─────────────────────────────────────────────────
#
# Each rule uses either:
#   container: "docker://quay.io/biocontainers/..."  (single bioinformatics tool)
#   conda:     "envs/<env>.yaml"                     (Python/R multi-package rules)
#
# Containers are pulled on first use via Apptainer (no docker build required).
# Conda environments are created on first use by Snakemake.
#
# Prerequisites:
#   - Apptainer (https://apptainer.org/docs/admin/main/installation.html)
#   - Conda / Mamba (https://conda-forge.org/miniforge/)
#
# Usage:
#   source setup.sh
#   integration_mapping_run WORKDIR/ 8
#
# The genome FASTA (00.genome/hg38.fna) must be BWA-indexed before the first run.
# Index it with:
#   integration_mapping_index_genome WORKDIR/
# ──────────────────────────────────────────────────────────────────────────────

integration_mapping_run_func() {
    WORKDIR=$(realpath ${1})
    THREADS=${2:-8}

    snakemake \
        --snakefile "${INTEGRATION_MAPPING_REPO_DIR}/snakemake/Snakefile" \
        --software-deployment-method apptainer conda \
        --apptainer-args "--bind ${WORKDIR}:/integration-mapping-pipeline/WORKDIR" \
        --directory "${INTEGRATION_MAPPING_REPO_DIR}/snakemake" \
        -j "${THREADS}" \
        --keep-going
}

integration_mapping_index_genome_func() {
    WORKDIR=$(realpath ${1})
    GENOME="${WORKDIR}/00.genome/hg38.fna"

    apptainer exec \
        --bind "${WORKDIR}:/integration-mapping-pipeline/WORKDIR" \
        "docker://quay.io/biocontainers/bwa:0.7.17--hed695b0_7" \
        bwa index "${GENOME}"
}

alias integration_mapping_run="integration_mapping_run_func"
alias integration_mapping_index_genome="integration_mapping_index_genome_func"
