# Integration Mapping Pipeline — Methods

This document describes the analysis performed by the Snakemake workflow in
`snakemake/`. Each section covers one or more pipeline steps and notes the
criteria that determine which reads, events, or sites are retained. Parameters
not mentioned are standard tool defaults.

---

## 1. UMI extraction and read cleaning (`clean_umi`)

Each R2 read contains a unique molecular identifier (UMI) embedded between the
sequencing primer and the attD sequence. Per-sample parameters (primer sequence,
UMI coordinates within the read, and optional i7 stagger) are drawn from
`metadata.tsv`.

- If the primer sequence is found in R2, the read is trimmed to begin at the
  primer. Otherwise, the stagger length is trimmed from the 5′ end.
- The 12 bp UMI is extracted from R2 at positions `[umi_start:umi_end]`
  (0-indexed start, 1-indexed end) and appended to both read headers as
  `:UMI_<sequence>` for downstream deduplication.
- The UMI region in R2 is replaced with poly-A, and its reverse complement in
  R1 is replaced with poly-T, to prevent the UMI sequence from interfering with
  alignment.

No reads are discarded at this step.

---

## 2. Adapter trimming (`fastp`)

Paired-end reads are trimmed with fastp using hardcoded adapter sequences:
R1 TruSeq (`AGATCGGAAGAGCACACGTCTGAACTCCAGTCA`) and R2 Nextera
(`CTGTCTCTTATACACATCTGACGCTGCCGACGA`). Standard quality trimming is applied.

---

## 3. Independent read alignments

Three independent alignments are performed on trimmed R1 and R2 reads
separately (not as pairs), each producing per-read BAM files:

- **Donor alignment** (`bwa mem`): each read is aligned to its
  sample-specific donor plasmid FASTA. Unmapped reads are discarded (`-F 4`).
- **Human alignment** (`bwa mem`): each read is aligned to the hg38 analysis
  set reference genome. Unmapped reads are discarded.
- **Donor-check alignment** (`bwa mem -a`): each read is aligned to a
  multi-sequence FASTA (`donor_check.fna`) containing a representative sequence
  for every donor construct in the experiment. All alignments are retained
  (`-a`), as this step is used for cross-contamination detection (see §4).

---

## 4. Per-read classification (`analyze_reads`)

This is the central classification step. Every trimmed read pair is assigned
values for several attributes that are used as gates in subsequent steps.

### 4.1 Donor verification (`r1_donor_check`, `r2_donor_check`)

For each read, the donor-check BAM is searched for the alignment with the
highest score (`percent_identity × alignment_length`). If ties occur, all
tied reference names are retained. The reference name(s) are split on `|` to
recover individual sample identifiers. A read is marked `True` if the
expected sample identifier appears among the top-scoring donors, and `False`
otherwise. Reads absent from the donor-check BAM are marked `UNKNOWN`.

This step detects cross-contamination: reads from sample A that fortuitously
align to sample B's donor are marked `False` for sample A and excluded
downstream.

**Note on naming:** the sample identifier used here is the full FASTQ filename
stem (e.g. `Em12-biorep2-techrep2`), which must exactly match the FASTA
sequence headers in `donor_check.fna`. Short LSR names (e.g. `Em12`) will
cause all reads to be marked `False` and silently empty the junction output.

### 4.2 Alignment category (`r1_align`, `r2_align`)

Each read is assigned one of four categories:

| Category | Meaning |
|---|---|
| `GENOMIC_ONLY` | Aligns to human genome, not to donor |
| `DONOR_ONLY` | Aligns to donor, not to human genome |
| `GENOMIC+DONOR` | Aligns to both |
| `NOTHING` | Unmapped to both |

### 4.3 Primer status (`primer_status`)

If R2 aligns to the donor, its alignment start position is compared to the
primer start coordinate from the donor BED file. If the alignment starts within
±2 bp of the primer position, the pair is marked `PRIMER`; otherwise
`NOPRIMER`. Pairs where R2 does not map to the donor are `NA`.

Only `PRIMER`-marked pairs are carried forward as candidate integration events.

### 4.4 attD integrity status (`r1_attd_status`, `r2_attd_status`)

If a read aligns to the donor, its coverage relative to the attD region
coordinates (from the donor BED) is classified as:

| Status | Meaning |
|---|---|
| `INTACT` | Read spans the full attD (refstart < attd_start and refend > attd_end) |
| `CLIPPED_LEFT` | Soft-clipped on the left within attD (clip > 4 bp) |
| `CLIPPED_RIGHT` | Soft-clipped on the right within attD (clip > 4 bp) |
| `UNKNOWN` | Other configurations |
| `NA` | Read does not map to donor |

Junction reads are expected to have clipped, not intact, attD alignments.
The 4 bp soft-clip threshold excludes marginal alignment artifacts.

### 4.5 Blacklist status (`r1_genome_blacklist`, `r2_genome_blacklist`)

The genomic alignment position of each read is intersected with the ENCODE
hg38 blacklist (`00.blacklist/hg38.bed`). Reads overlapping blacklisted
regions carry the region name; others carry `None` or `NA`.

---

## 5. Junction read extraction (`get_junction_reads`)

Reads are selected as junction reads — spanning a donor–genome integration
boundary — only if **all** of the following conditions are satisfied:

| Condition | Rationale |
|---|---|
| `r2_donor_check == "True"` | R2 confirms alignment to the correct donor (not contamination) |
| `r1_donor_check != "True"` | R1 is not donor-confirmed (it should be the genomic arm) |
| `r1_align` contains `"GENOMIC"` | R1 has a genomic alignment component |
| `primer_status == "PRIMER"` | Integration-associated primer is present |
| Neither attD status is `INTACT` | attD is clipped, consistent with a junction |
| Not both reads `GENOMIC_ONLY` | At least one read contacts the donor |
| Not both reads `DONOR_ONLY` | At least one read contacts the genome |
| No `ONLY` + `NOTHING` pairing | Both reads must have an alignment |
| `r1_genome_blacklist` is `"None"` or `"NA"` | R1 not in a blacklisted region |
| `r2_genome_blacklist` is `"None"` or `"NA"` | R2 not in a blacklisted region |

Reads failing any criterion are excluded. The result is a FASTQ of chimeric
read pairs whose R1 spans the integration junction in the host genome.

---

## 6. Donor read extraction (`get_donor_reads`)

A parallel branch selects read pairs that map entirely within the donor
construct (no genomic component). These are used to quantify donor plasmid
abundance independently of integration events.

Criteria: `r2_donor_check == "True"`, `r1_donor_check != "True"`,
`primer_status == "PRIMER"`, both reads `DONOR_ONLY`, and neither attD status
is `CLIPPED_LEFT` or `CLIPPED_RIGHT` (intact attD required — the opposite
criterion from junction reads).

---

## 7. Paired re-alignment of junction reads to human genome (`bwa_align_both_human`)

Junction read pairs are re-aligned together (`bwa mem`) to hg38, then filtered:

- **MAPQ ≥ 30**: retains only high-confidence, uniquely mapping alignments.
- **Supplementary alignments removed** (`-F 2048`).
- **Insert size filter**: read pairs with |TLEN| ≥ 1500 bp are discarded as
  likely mapping artifacts.
- **Blacklist removal**: aligned reads are further filtered with
  `bedtools intersect -v` against the blacklist BED.

PCR duplicates are then removed with `samtools markdup -r` after collating,
fixing mate information, and coordinate-sorting (the markdup BAM is used for
unique read counts).

---

## 8. Integration site detection (`mgefinder`)

The deduplicated junction BAM is passed to MGEfinder (`mgefinder find`) with
`-mcc 1 -mincount 1`, which detects candidate integration positions from
soft-clipped read clusters with no minimum count threshold. All positions
supported by at least one soft-clipped read are reported.

---

## 9. Breakpoint resolution (`kmer_flanks`)

Each MGEfinder candidate is refined to a precise breakpoint using two
independent methods applied to the soft-clipped consensus sequence:

**K-mer method:** Every 9-mer in the consensus is searched against the attD
sequence. The genomic cut position is inferred from the most frequently
occurring match position across all 9-mers.

**Pairwise alignment method:** A local Smith–Waterman alignment (match=2,
mismatch=−1, gap-open=−2, gap-extend=−2, no end-gap penalties) is performed
between the consensus and the attD sequence. The alignment boundaries are
mapped back to genomic coordinates.

Both methods produce an `align_cut_dist` (distance from the predicted
breakpoint to the soft-clip position) and alignment identity metrics
(`pident_trimmed`, `pident_contig`). Sites are assigned to loci by merging
overlapping ±500 bp windows around each candidate position.

---

## 10. Read and UMI counting per locus (`raw_counts`, `umi_counts`)

For each locus, reads in the junction BAM are counted in two ways:

- **RAW**: all reads in the unfiltered `human.both.bam` that intersect the
  locus (via `bedtools intersect`).
- **UNIQ**: reads from the duplicate-marked BAM (`human.both.markdup.bam`),
  representing PCR-deduplicated read count.

UMI-based deduplication is performed separately: reads intersecting each
locus have their UMI extracted from the read header. UMIs are filtered to
exactly 12 bp and at most 1 ambiguous (`N`) base, then clustered using the
directional method (edit-distance threshold = 1) to collapse PCR duplicates
sharing similar UMIs. The cluster count is the UMI count for that locus.

For donor reads, a parallel UMI count is computed from the
`bwa_align_both_donor` BAM. Only forward-strand reads from proper pairs with
fragment length ≥ 100 bp are included to avoid double-counting.

---

## 11. Integration site curation (`integration_sites`)

Candidate loci are filtered to produce the final set of high-confidence
integration sites:

| Filter | Threshold | Rationale |
|---|---|---|
| Breakpoint precision | −15 ≤ `align_cut_dist` ≤ 15 bp | Excludes imprecisely localized sites |
| Alignment identity | `pident_trimmed` > 0.80 **or** `pident_contig` > 0.80 | Requires ≥80% sequence match to attD |

Passing sites are ranked by UNIQ read count (descending), then soft-clip
count (descending). Within each locus, only the top-ranked site is retained
(deduplication). A maximum of 200 sites per sample are carried into the
final output.

---

## 12. Genomic annotation (`get_site_annotations`)

Each curated site is intersected with GENCODE v38 gene and exon annotations
(filtered to `gene` and `exon` features respectively, sorted to match the
chromosome order of `hg38.genome.txt`) using `bedtools intersect -sorted`.
Sites are classified as:

- **EXONIC**: overlaps a gene and at least one exon of that gene.
- **INTRONIC**: overlaps a gene but no exon.
- **INTERGENIC**: no gene overlap; chromosome name ≤ 5 characters (canonical chromosomes).
- **UNANNOTATED**: no gene overlap on an unplaced or alternate contig (chromosome name > 5 characters).

Overlapping gene IDs are recorded for each site.

---

## 13. Finalization (`finalize_counts`)

Per-sample count tables, site sequences, and annotations are merged across
all samples using `finalize_counts.R`. Integration site sequences (±30 bp
flanking) are consolidated into `all_integration_sites.fna` by
`extract_integration_sites.py`. Sequence logos for each sample's integration
sites are generated with `seqlogos.R`.

---

## Summary of filtering thresholds

| Stage | Parameter | Value |
|---|---|---|
| Primer position tolerance | ±2 bp | `analyze_reads` |
| attD soft-clip threshold | > 4 bp | `analyze_reads` |
| Junction re-alignment MAPQ | ≥ 30 | `bwa_align_both_human/donor` |
| Junction insert size | \|TLEN\| < 1500 bp | `bwa_align_both_human` |
| MGEfinder minimum read count | 1 | `mgefinder` |
| K-mer size for breakpoint mapping | 9 bp | `kmer_flanks` |
| Breakpoint precision filter | ±15 bp | `integration_sites` |
| Sequence identity filter | > 80% | `integration_sites` |
| Maximum sites per sample | 200 | `integration_sites` |
| UMI length requirement | exactly 12 bp | `umi_counts`, `donor_umi_counts` |
| UMI ambiguous base tolerance | ≤ 1 N | `umi_counts` |
| UMI clustering edit distance | 1 | `umi_counts`, `donor_umi_counts` |
| Donor read minimum fragment length | ≥ 100 bp | `donor_umi_counts` |
