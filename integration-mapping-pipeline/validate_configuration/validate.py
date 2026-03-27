#!/usr/bin/env python3
"""
Validates an integration-mapping-pipeline workdir configuration.

Exit codes:
  0 — all checks passed (VALID)
  1 — one or more checks failed (INVALID)

Usage:
  python3 validate.py <workdir>
"""

import sys
import os
import csv
from pathlib import Path

errors = []
warnings = []


def err(msg):
    errors.append(msg)


def warn(msg):
    warnings.append(msg)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def check_file(path, label):
    """Return True and print nothing if path exists; record error and return False otherwise."""
    if not path.exists():
        err(f"Missing {label}: {path}")
        return False
    return True


def fai_chromosomes(fai_path):
    """Return ordered list of chromosome names from a .fai file."""
    chroms = []
    with open(fai_path) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if parts:
                chroms.append(parts[0])
    return chroms


def genome_txt_chromosomes(genome_txt_path):
    """Return ordered list of chromosome names from a genome.txt (two-column TSV)."""
    chroms = []
    with open(genome_txt_path) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if parts:
                chroms.append(parts[0])
    return chroms


def gff3_chromosomes(gff3_path):
    """Return ordered list of chromosome names (deduplicated, in first-occurrence order) from a GFF3."""
    seen = []
    seen_set = set()
    with open(gff3_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            chrom = line.split("\t")[0]
            if chrom not in seen_set:
                seen.append(chrom)
                seen_set.add(chrom)
    return seen


def fasta_headers(fasta_path):
    """Return list of full header strings (without '>') from a FASTA file."""
    headers = []
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                headers.append(line[1:].strip())
    return headers


def bed_feature_names(bed_path):
    """Return set of feature names from column 4 of a BED file."""
    names = set()
    with open(bed_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 4:
                names.add(parts[3])
    return names


# ---------------------------------------------------------------------------
# Check 1 — Workdir exists
# ---------------------------------------------------------------------------

def check_workdir_exists(wd):
    if not wd.exists():
        err(f"Workdir does not exist: {wd}")
        return False
    if not wd.is_dir():
        err(f"Workdir path is not a directory: {wd}")
        return False
    return True


# ---------------------------------------------------------------------------
# Check 2 — Required directories
# ---------------------------------------------------------------------------

REQUIRED_DIRS = [
    "00.blacklist",
    "00.donor_check",
    "00.donor_map",
    "00.fastq",
    "00.genome",
]


def check_required_dirs(wd):
    ok = True
    for d in REQUIRED_DIRS:
        path = wd / d
        if not path.is_dir():
            err(f"Missing required directory: {path}")
            ok = False
    return ok


# ---------------------------------------------------------------------------
# Check 3 — Required files
# ---------------------------------------------------------------------------

def check_required_files(wd):
    ok = True
    required = {
        "metadata.tsv": wd / "metadata.tsv",
        "hg38 blacklist BED": wd / "00.blacklist" / "hg38.bed",
        "donor_check.fna": wd / "00.donor_check" / "donor_check.fna",
        "hg38.fna": wd / "00.genome" / "hg38.fna",
        "hg38.fna.fai": wd / "00.genome" / "hg38.fna.fai",
        "hg38.genome.txt": wd / "00.genome" / "hg38.genome.txt",
        "hg38.genes.gff3": wd / "00.genome" / "hg38.genes.gff3",
        "hg38.exons.gff3": wd / "00.genome" / "hg38.exons.gff3",
    }
    for label, path in required.items():
        if not check_file(path, label):
            ok = False

    # BWA index files
    for ext in ("amb", "ann", "bwt", "pac", "sa"):
        p = wd / "00.genome" / f"hg38.fna.{ext}"
        if not p.exists():
            err(f"Missing BWA index file: {p}")
            ok = False

    return ok


# ---------------------------------------------------------------------------
# Check 4 — metadata.tsv columns
# ---------------------------------------------------------------------------

REQUIRED_METADATA_COLS = {
    "sample", "lsr", "primer_seq", "umi_start", "umi_end",
}


def check_metadata(wd):
    meta_path = wd / "metadata.tsv"
    if not meta_path.exists():
        return None  # already reported in check 3

    samples = []
    try:
        with open(meta_path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            cols = set(reader.fieldnames or [])
            missing_cols = REQUIRED_METADATA_COLS - cols
            if missing_cols:
                err(f"metadata.tsv missing required columns: {sorted(missing_cols)}")
                return None
            for row in reader:
                s = row.get("sample", "").strip()
                if s:
                    samples.append(s)
    except Exception as e:
        err(f"Could not parse metadata.tsv: {e}")
        return None

    if not samples:
        err("metadata.tsv contains no sample rows")
        return None

    return samples


# ---------------------------------------------------------------------------
# Check 5 — FASTQs present and paired; no orphans
# ---------------------------------------------------------------------------

def check_fastqs(wd, metadata_samples):
    fastq_dir = wd / "00.fastq"
    if not fastq_dir.is_dir():
        return  # already reported

    r1_files = {p.name.replace(".R1.fq.gz", "") for p in fastq_dir.glob("*.R1.fq.gz")}
    r2_files = {p.name.replace(".R2.fq.gz", "") for p in fastq_dir.glob("*.R2.fq.gz")}

    for sample in (metadata_samples or []):
        if sample not in r1_files:
            err(f"Missing R1 FASTQ for sample '{sample}': expected {fastq_dir}/{sample}.R1.fq.gz")
        if sample not in r2_files:
            err(f"Missing R2 FASTQ for sample '{sample}': expected {fastq_dir}/{sample}.R2.fq.gz")

    # Orphaned R1 without R2 and vice-versa
    for s in r1_files - r2_files:
        err(f"R1 FASTQ has no matching R2: {s}.R1.fq.gz")
    for s in r2_files - r1_files:
        err(f"R2 FASTQ has no matching R1: {s}.R2.fq.gz")

    # FASTQs in directory with no metadata entry
    if metadata_samples is not None:
        meta_set = set(metadata_samples)
        for s in r1_files - meta_set:
            warn(f"FASTQ '{s}' has no entry in metadata.tsv")


# ---------------------------------------------------------------------------
# Check 6 — donor_map files per sample
# ---------------------------------------------------------------------------

EXPECTED_BED_FEATURES = {"Primer", "UMI", "attD"}


def check_donor_map(wd, metadata_samples):
    donor_map = wd / "00.donor_map"
    if not donor_map.is_dir():
        return

    for sample in (metadata_samples or []):
        fasta = donor_map / f"{sample}.fasta"
        bed = donor_map / f"{sample}.bed"
        attd = donor_map / f"{sample}.attd.fasta"

        for p, label in [(fasta, "FASTA"), (bed, "BED"), (attd, "attD FASTA")]:
            if not p.exists():
                err(f"Missing donor_map {label} for sample '{sample}': {p}")

        if bed.exists():
            features = bed_feature_names(bed)
            missing = EXPECTED_BED_FEATURES - features
            if missing:
                err(
                    f"donor_map BED for '{sample}' is missing feature(s) {sorted(missing)}; "
                    f"found: {sorted(features)}"
                )


# ---------------------------------------------------------------------------
# Check 7 — donor_check.fna headers cover every sample exactly once
# ---------------------------------------------------------------------------

def check_donor_check(wd, metadata_samples):
    fna = wd / "00.donor_check" / "donor_check.fna"
    if not fna.exists():
        return  # already reported

    headers = fasta_headers(fna)

    # Build mapping: sample_name -> list of header strings that claim it
    sample_to_headers = {}
    for h in headers:
        # Each header may be pipe-separated: Name1|Name2|...
        for part in h.split("|"):
            part = part.strip()
            sample_to_headers.setdefault(part, []).append(h)

    if metadata_samples is None:
        return

    for sample in metadata_samples:
        if sample not in sample_to_headers:
            err(
                f"Sample '{sample}' not found in any donor_check.fna header. "
                f"Headers must use the full FASTQ filename stem (e.g. '{sample}'), "
                f"not a short LSR name. This will cause all reads to be marked "
                f"r1/r2_donor_check=False and silently empty junction output."
            )
        elif len(sample_to_headers[sample]) > 1:
            err(
                f"Sample '{sample}' appears in multiple donor_check.fna headers: "
                f"{sample_to_headers[sample]}"
            )


# ---------------------------------------------------------------------------
# Check 8 — hg38.genome.txt chromosome order matches hg38.fna.fai
# ---------------------------------------------------------------------------

def check_genome_txt_order(wd):
    fai = wd / "00.genome" / "hg38.fna.fai"
    genome_txt = wd / "00.genome" / "hg38.genome.txt"
    if not fai.exists() or not genome_txt.exists():
        return  # already reported

    fai_chroms = fai_chromosomes(fai)
    txt_chroms = genome_txt_chromosomes(genome_txt)

    if fai_chroms != txt_chroms:
        # Find first mismatch for a helpful message
        for i, (a, b) in enumerate(zip(fai_chroms, txt_chroms)):
            if a != b:
                err(
                    f"hg38.genome.txt chromosome order does not match hg38.fna.fai. "
                    f"First mismatch at position {i}: fai has '{a}', genome.txt has '{b}'. "
                    f"Regenerate genome.txt with: cut -f1,2 hg38.fna.fai > hg38.genome.txt"
                )
                return
        # One is longer than the other
        if len(fai_chroms) != len(txt_chroms):
            err(
                f"hg38.genome.txt has {len(txt_chroms)} chromosomes but hg38.fna.fai has "
                f"{len(fai_chroms)}. Regenerate with: cut -f1,2 hg38.fna.fai > hg38.genome.txt"
            )


# ---------------------------------------------------------------------------
# Check 9 — GFF3 files sorted to match genome.txt
# ---------------------------------------------------------------------------

def check_gff3_order(wd):
    genome_txt = wd / "00.genome" / "hg38.genome.txt"
    if not genome_txt.exists():
        return  # already reported

    txt_chroms = genome_txt_chromosomes(genome_txt)
    txt_order = {c: i for i, c in enumerate(txt_chroms)}

    for gff_name in ("hg38.genes.gff3", "hg38.exons.gff3"):
        gff = wd / "00.genome" / gff_name
        if not gff.exists():
            continue  # already reported

        gff_chroms = gff3_chromosomes(gff)
        # Filter to only chromosomes that appear in genome.txt
        known = [c for c in gff_chroms if c in txt_order]
        indices = [txt_order[c] for c in known]

        # Check that indices are non-decreasing
        for i in range(1, len(indices)):
            if indices[i] < indices[i - 1]:
                err(
                    f"{gff_name} is not sorted to match hg38.genome.txt. "
                    f"Chromosome '{known[i]}' (genome.txt position {indices[i]}) appears after "
                    f"'{known[i-1]}' (position {indices[i-1]}). "
                    f"Re-sort with: bedtools sort -g hg38.genome.txt -i {gff_name} > sorted.gff3 "
                    f"&& mv sorted.gff3 {gff_name}"
                )
                break


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) < 2:
        print("Usage: validate.py <workdir>", file=sys.stderr)
        sys.exit(2)

    wd = Path(sys.argv[1]).resolve()

    if not check_workdir_exists(wd):
        print("INVALID")
        for e in errors:
            print(f"  ERROR: {e}")
        sys.exit(1)

    check_required_dirs(wd)
    check_required_files(wd)
    metadata_samples = check_metadata(wd)
    check_fastqs(wd, metadata_samples)
    check_donor_map(wd, metadata_samples)
    check_donor_check(wd, metadata_samples)
    check_genome_txt_order(wd)
    check_gff3_order(wd)

    if warnings:
        for w in warnings:
            print(f"  WARNING: {w}")

    if errors:
        print("INVALID")
        for e in errors:
            print(f"  ERROR: {e}")
        sys.exit(1)
    else:
        print("VALID")
        sys.exit(0)


if __name__ == "__main__":
    main()
