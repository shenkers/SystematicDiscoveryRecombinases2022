# Reconstituting the 3-sample integration mapping test

This folder contains everything needed to reconstitute a 3-sample subset of
the integration mapping pipeline on a fresh EC2 instance. Large files (FASTQs,
genome) are retrieved in the steps below; all other required files are already
present here.

All shell snippets below use `PROJECT_ROOT` to refer to the directory where you
want to run the pipeline. Set it once before following the steps:

```bash
export PROJECT_ROOT="/path/to/your/working/directory"
```

## Samples

| Pipeline sample name      | SRX accession | SRR accession |
|---------------------------|---------------|---------------|
| Cp36-biorep1-techrep1     | SRX17313194   | SRR21306543   |
| Em12-biorep2-techrep2     | SRX17313211   | SRR21306526   |
| Enc9-biorep1-techrep1     | SRX17313212   | SRR21306525   |

All three are from BioProject **PRJNA778877**.

---

## Step 1 — Download FASTQs with nf-core/fetchngs

`ids.csv` in this folder lists the three SRR accessions. Nextflow is assumed
to already be installed.

```bash
cd /path/to/workdir          # wherever you want to run from
nextflow run nf-core/fetchngs \
    --input ${PROJECT_ROOT}/ids.csv \
    --outdir ./fetchngs_out \
    -profile docker \
    -r 1.12.0
```

fetchngs will deposit paired FASTQs named `{SRX}_{SRR}_1.fastq.gz` /
`{SRX}_{SRR}_2.fastq.gz` under `./fetchngs_out/fastq/`.

---

## Step 2 — Rename/link FASTQs into the workdir

The pipeline expects files named `{sample}.R1.fq.gz` / `{sample}.R2.fq.gz`
in `workdir/00.fastq/`. Create symlinks with the right names:

```bash
FETCHNGS_FASTQ="./fetchngs_out/fastq"
FASTQ_DIR="${PROJECT_ROOT}/workdir/00.fastq"

ln -sf "${FETCHNGS_FASTQ}/SRX17313194_SRR21306543_1.fastq.gz" "${FASTQ_DIR}/Cp36-biorep1-techrep1.R1.fq.gz"
ln -sf "${FETCHNGS_FASTQ}/SRX17313194_SRR21306543_2.fastq.gz" "${FASTQ_DIR}/Cp36-biorep1-techrep1.R2.fq.gz"

ln -sf "${FETCHNGS_FASTQ}/SRX17313211_SRR21306526_1.fastq.gz" "${FASTQ_DIR}/Em12-biorep2-techrep2.R1.fq.gz"
ln -sf "${FETCHNGS_FASTQ}/SRX17313211_SRR21306526_2.fastq.gz" "${FASTQ_DIR}/Em12-biorep2-techrep2.R2.fq.gz"

ln -sf "${FETCHNGS_FASTQ}/SRX17313212_SRR21306525_1.fastq.gz" "${FASTQ_DIR}/Enc9-biorep1-techrep1.R1.fq.gz"
ln -sf "${FETCHNGS_FASTQ}/SRX17313212_SRR21306525_2.fastq.gz" "${FASTQ_DIR}/Enc9-biorep1-techrep1.R2.fq.gz"
```

---

## Step 3 — Download and prepare the hg38 reference genome

The pipeline uses the **UCSC hg38 analysis set** (hard-masked, `chr`-prefixed
chromosome names). Download and concatenate the per-chromosome FASTAs, then
build the BWA index and supporting files.

```bash
GENOME_DIR="${PROJECT_ROOT}/workdir/00.genome"
cd "${GENOME_DIR}"

# Download the per-chromosome FASTA archive (~900 MB)
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.chroms.tar.gz

# Extract and concatenate into a single FASTA
tar -xzf hg38.analysisSet.chroms.tar.gz
cat hg38.analysisSet.chroms/chr*.fa > hg38.fna
rm -rf hg38.analysisSet.chroms hg38.analysisSet.chroms.tar.gz

# Build BWA index (~1 hour, ~9 GB)
bwa index hg38.fna

# FASTA index and chromosome sizes
samtools faidx hg38.fna
cut -f1,2 hg38.fna.fai > hg38.genome.txt
```

---

## Step 4 — Download and filter GENCODE v38 gene annotations

The pipeline needs two GFF3 files filtered from the GENCODE v38 annotation:
one containing only `gene` features and one containing only `exon` features.

**Important:** The pipeline uses `bedtools intersect -sorted` with these GFF3 files
and with the aligned BAM files, which requires all of them to share exactly the
same chromosome order as `hg38.genome.txt`. Two constraints apply:

1. `hg38.genome.txt` must be derived directly from `hg38.fna.fai` (see Step 3)
   so that it matches the chromosome order of the FASTA and the BAMs that were
   aligned to it. Do not substitute a pre-downloaded chromosome sizes file — they
   often use a different ordering (e.g. sorted by size) that will not match the
   BAM headers.

2. GENCODE GFF3 files use lexicographic chromosome order (chr1, chr10, chr11, …
   chr2, …) while the hg38 analysis set FASTA, when concatenated with
   `cat chr*.fa`, also uses lexicographic order — these happen to match, but the
   sort step below is still required if you obtained genome.txt from any other
   source. Skipping it will cause `bedtools` to abort with a "different sort
   order" error at the `get_site_annotations`, `raw_counts`, or `umi_counts`
   steps.

```bash
GENOME_DIR="${PROJECT_ROOT}/workdir/00.genome"
cd "${GENOME_DIR}"

# Download GENCODE v38 annotation (~50 MB compressed)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gff3.gz

# Filter gene and exon features (preserve header lines beginning with #)
zcat gencode.v38.annotation.gff3.gz | awk '$3=="gene" || /^#/' > hg38.genes.gff3
zcat gencode.v38.annotation.gff3.gz | awk '$3=="exon" || /^#/' > hg38.exons.gff3

rm gencode.v38.annotation.gff3.gz

# Sort GFF3 files to match the chromosome order in hg38.genome.txt.
# bedtools sort preserves lines beginning with '#' at the top.
bedtools sort -g hg38.genome.txt -i hg38.genes.gff3 > hg38.genes.sorted.gff3 \
    && mv hg38.genes.sorted.gff3 hg38.genes.gff3
bedtools sort -g hg38.genome.txt -i hg38.exons.gff3 > hg38.exons.sorted.gff3 \
    && mv hg38.exons.sorted.gff3 hg38.exons.gff3
```

---

## Critical: donor_check.fna sequence headers must match FASTQ sample names

The pipeline determines which donor plasmid to expect for each sample using
`params.sample`, which is derived directly from the FASTQ filename stem — i.e.
the part before `.R1.fq.gz`. For example, a file named
`Em12-biorep2-techrep2.R1.fq.gz` produces `params.sample = Em12-biorep2-techrep2`.

The pipeline verifies each read by aligning it to `00.donor_check/donor_check.fna`
and checking that the best-matching reference name equals `params.sample`. This
means the FASTA headers in `donor_check.fna` **must exactly match the FASTQ
filename stems** — not the short LSR name from `metadata.tsv`.

If multiple samples in the same run share an identical donor sequence (same LSR,
different biorep/techrep), they should be combined into a single FASTA record
with names separated by `|`:

```
# correct — headers match the FASTQ stems
>Em12-biorep1-techrep1|Em12-biorep2-techrep2
GACGGGCACC...

# wrong — short LSR name never matches params.sample
>Em12
GACGGGCACC...
```

A mismatch here causes `r1_donor_check` and `r2_donor_check` to be `False` for
every read, which silently empties the junction read output (step 05) and all
downstream steps, with no obvious error message.

---

## Step 5 — Verify the workdir layout

After completing the steps above, `workdir/` should look like this:

```
workdir/
├── metadata.tsv                        ← included
├── 00.blacklist/
│   └── hg38.bed                        ← included
├── 00.donor_check/
│   └── donor_check.fna                 ← included
├── 00.donor_map/
│   ├── Cp36-biorep1-techrep1.fasta     ← included
│   ├── Cp36-biorep1-techrep1.bed       ← included
│   ├── Cp36-biorep1-techrep1.attd.fasta← included
│   ├── Em12-biorep2-techrep2.*         ← included
│   └── Enc9-biorep1-techrep1.*         ← included
├── 00.fastq/
│   ├── Cp36-biorep1-techrep1.R1.fq.gz  ← symlink from fetchngs (Step 2)
│   ├── Cp36-biorep1-techrep1.R2.fq.gz
│   ├── Em12-biorep2-techrep2.R{1,2}.fq.gz
│   └── Enc9-biorep1-techrep1.R{1,2}.fq.gz
└── 00.genome/
    ├── hg38.fna                         ← downloaded (Step 3)
    ├── hg38.fna.{amb,ann,bwt,pac,sa}    ← BWA index (Step 3)
    ├── hg38.fna.fai                     ← samtools index (Step 3)
    ├── hg38.genome.txt                  ← derived (Step 3)
    ├── hg38.genes.gff3                  ← derived (Step 4)
    └── hg38.exons.gff3                  ← derived (Step 4)
```

---

## Step 6 — Run the pipeline

Build the Docker image (once) from the repo:

```bash
cd /path/to/SystematicDiscoveryRecombinases2022/integration-mapping-pipeline
docker build -t integration_mapping env/
```

Then run:

```bash
docker run --rm \
    -v /path/to/SystematicDiscoveryRecombinases2022/integration-mapping-pipeline/snakemake:/integration-mapping-pipeline/snakemake \
    -v ${PROJECT_ROOT}/workdir:/integration-mapping-pipeline/WORKDIR \
    integration_mapping \
    snakemake -j 8 --keep-going --config wd=/integration-mapping-pipeline/WORKDIR
```

Results will appear under `workdir/11.results/`. See the pipeline README for a
description of output files.
