#!/usr/bin/env Rscript
# =====================================================
# sequence_utils.R
#
# Small header / FASTA utilities shared by the alignment and QC helpers:
#   - get_accession()       parse a GenBank accession out of a FASTA header
#   - standardize_taxon()   parse a header down to a "Genus_species" label
#   - bad_fraction()        fraction of gap/ambiguous characters in a sequence
#   - read_combined_locus() read + merge the ingroup and outgroup FASTA for a locus
#   - subset_to_species()   keep only sequences whose binomial is in a given set
#
# Inputs : per-locus FASTA files under <INGROUP_DIR>/<locus>/ and <OUTGROUP_DIR>/<locus>/.
# Outputs: none (pure helpers; operate on in-memory DNAStringSets).
#
# Provenance:
#   Written by N. Weksler for the consolidated Danaini tree-building workflow. These consolidate
#   and fix the header/accession parsing that was duplicated (with three inconsistent accession
#   regexes) across Marianna Karageorgi's original per-step scripts.
# =====================================================

suppressPackageStartupMessages({
  library(Biostrings)  # DNAStringSet I/O
  library(stringr)     # str_extract
})

# Unified accession regex — matches RefSeq (XM_032672544.2) AND standard (EU141321.1);
# falls back to the first header token when nothing matches.
get_accession <- function(hdr) {
  acc <- str_extract(hdr, "[A-Z]{1,3}_?\\d{5,}(?:\\.\\d+)?")
  ifelse(is.na(acc), sub("^(\\S+).*$", "\\1", hdr), acc)
}

# Header -> Genus_species (NA if no binomial can be parsed).
standardize_taxon <- function(hdr) {
  sp <- str_extract(hdr, "[A-Z][a-z]+\\s+[a-z]+")
  if (is.na(sp)) return(NA_character_)
  gsub("\\s+", "_", sp)
}

# Fraction of gap / ambiguous characters in a single sequence.
bad_fraction <- function(x) {
  s <- strsplit(as.character(x), "", fixed = TRUE)[[1]]
  mean(s %in% c("-", "N", "n", "?"))
}

# Read the ingroup + outgroup FASTA for a locus and merge them in memory.
# INGROUP_DIR / OUTGROUP_DIR are expected to exist in the calling environment.
read_combined_locus <- function(locus, ingroup_dir = INGROUP_DIR, outgroup_dir = OUTGROUP_DIR) {
  f_in  <- file.path(ingroup_dir,  locus, paste0(locus, "_all_trimmed.fasta"))
  f_out <- file.path(outgroup_dir, locus, paste0(locus, "_all_trimmed.fasta"))
  seqs <- DNAStringSet()
  if (file.exists(f_in))  seqs <- c(seqs, readDNAStringSet(f_in))
  if (file.exists(f_out)) seqs <- c(seqs, readDNAStringSet(f_out))
  seqs
}

# Keep only sequences whose parsed binomial is in `keep_std` (a Genus_species vector).
subset_to_species <- function(dna, keep_std) {
  sp <- vapply(names(dna), standardize_taxon, character(1))
  dna[!is.na(sp) & sp %in% keep_std]
}
