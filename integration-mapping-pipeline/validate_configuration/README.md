# validate_configuration

Validates a workdir before running the integration-mapping-pipeline.

## Usage

```bash
bash validate.sh <path/to/workdir>
```

Prints a single line — `VALID` or `INVALID` — and exits 0 or 1 respectively.
When invalid, each violated expectation is printed as an `ERROR:` line before
the verdict. Potential issues that are not blocking are printed as `WARNING:`
lines (exit code is still 0 if there are no errors).

The script can be run from any directory; it locates `validate.py` relative
to itself and requires only Python 3 (no third-party packages).

## Example output

```
  ERROR: Sample 'Em12-biorep2-techrep2' not found in any donor_check.fna header. ...
  ERROR: hg38.genome.txt chromosome order does not match hg38.fna.fai. ...
INVALID
```

```
VALID
```

## What is checked

| # | Check | Why it matters |
|---|-------|----------------|
| 1 | Workdir path exists and is a directory | Fails immediately with a clear message instead of a cryptic downstream error |
| 2 | Required sub-directories present (`00.blacklist`, `00.donor_check`, `00.donor_map`, `00.fastq`, `00.genome`) | Pipeline rules expect these exact paths |
| 3 | Required files present (metadata.tsv, hg38.bed, donor_check.fna, hg38.fna + BWA index files, hg38.fna.fai, hg38.genome.txt, hg38.genes.gff3, hg38.exons.gff3) | Missing any of these causes a pipeline rule to fail |
| 4 | `metadata.tsv` has required columns and at least one sample row | Columns `sample`, `lsr`, `primer_seq`, `umi_start`, `umi_end` are read by `analyze_reads` |
| 5 | Every metadata sample has a paired R1 + R2 FASTQ; no orphaned FASTQs | Snakemake derives `SAMPLES` from the FASTQ glob; an orphan will produce a dead-end rule |
| 6 | Every metadata sample has `.fasta`, `.bed`, and `.attd.fasta` in `00.donor_map`; BED contains `Primer`, `UMI`, and `attD` features | These files drive donor alignment and attD/primer coordinate lookups in `analyze_reads` |
| 7 | Every metadata sample appears in exactly one `donor_check.fna` header using its **full** FASTQ filename stem | `params.sample` is the FASTQ stem (e.g. `Em12-biorep2-techrep2`); a short LSR name (e.g. `Em12`) or missing entry causes `r1/r2_donor_check` to be `False` for every read, silently emptying junction output with no error message |
| 8 | `hg38.genome.txt` chromosome order matches `hg38.fna.fai` | BAMs are aligned to hg38.fna, so they share its chromosome order; `bedtools intersect -sorted` requires genome.txt to match. Fix: `cut -f1,2 hg38.fna.fai > hg38.genome.txt` |
| 9 | `hg38.genes.gff3` and `hg38.exons.gff3` are sorted in the same chromosome order as `hg38.genome.txt` | `bedtools intersect -sorted` aborts with a "different sort order" error at `get_site_annotations`, `raw_counts`, or `umi_counts`. Fix: `bedtools sort -g hg38.genome.txt -i <file> > sorted && mv sorted <file>` |

## donor_check.fna header format

Headers in `donor_check.fna` must use the full FASTQ filename stem, not the
short LSR name from `metadata.tsv`. Multiple samples that share an identical
donor sequence can be combined with `|`:

```
# correct — header matches the FASTQ stems
>Em12-biorep1-techrep1|Em12-biorep2-techrep2
GACGGGCACC...

# wrong — short LSR name never matches params.sample
>Em12
GACGGGCACC...
```
