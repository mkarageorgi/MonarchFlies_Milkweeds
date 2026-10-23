#!/usr/bin/env Rscript
# =====================================================
# codon_alignment.R
#
# Codon-aware alignment of a protein-coding locus:
#   1) pick the best reading frame per sequence (minimize internal stops, maximize length)
#   2) translate to amino acids
#   3) align the amino acids with MAFFT (more reliable than nucleotide alignment)
#   4) back-translate the AA alignment to a codon-aligned nucleotide alignment
#
# Main entry point:
#   codon_align(dna, gcode, mafft_bin, threads) -> codon-aligned DNAStringSet
#
# Genetic codes:
#   - nuclear loci : Biostrings::GENETIC_CODE (standard)
#   - COI          : Biostrings::getGeneticCode("5")  (invertebrate mitochondrial)
#
# Provenance:
#   The codon-aware alignment logic (best-frame -> translate -> MAFFT on AAs -> back-translate)
#   is Marianna Karageorgi's, originally in codon_align_locus.R. N. Weksler refactored it to
#   operate entirely in memory on a DNAStringSet (rather than reading/writing per-locus FASTA
#   files), and to drop sequences with empty translations before MAFFT so they cannot corrupt
#   the alignment.
#
# =====================================================

suppressPackageStartupMessages({
  library(Biostrings)  # translate(), DNA/AA string sets
  library(stringr)     # str_count, str_sub
})

# --- internals (prefixed with "." ; not meant to be called directly) ---

.trim_to_full_codons <- function(s) {
  s <- toupper(gsub("\\s+", "", s)); n3 <- nchar(s) - (nchar(s) %% 3)
  if (n3 < 3) "" else substr(s, 1, n3)
}

.translate_frame <- function(dna_str, frame, gcode) {
  dna_str <- toupper(gsub("\\s+", "", dna_str))
  if (nchar(dna_str) < frame) return("")
  s <- .trim_to_full_codons(substr(dna_str, frame, nchar(dna_str)))
  if (nchar(s) < 3) return("")
  as.character(suppressWarnings(Biostrings::translate(DNAString(s),
      genetic.code = gcode, if.fuzzy.codon = "X")))
}

# Score a translation: fewer internal stops is better, then longer is better.
.score_tr <- function(aa) {
  if (nchar(aa) == 0) return(c(stops = Inf, len = 0))
  stops <- str_count(aa, fixed("*")); if (endsWith(aa, "*")) stops <- stops - 1
  c(stops = stops, len = nchar(aa))
}

# Choose the best of the three reading frames for one sequence.
.best_frame <- function(dna_str, gcode) {
  aa <- lapply(1:3, function(f) .translate_frame(dna_str, f, gcode))
  sc <- t(vapply(aa, .score_tr, numeric(2)))
  best <- order(sc[, "stops"], -sc[, "len"])[1]
  dna_in <- .trim_to_full_codons(substr(toupper(dna_str), best, nchar(dna_str)))
  list(aa = aa[[best]], dna = dna_in)
}

# Back-translate an AA alignment row to codons using the in-frame nucleotide sequence.
.backtranslate <- function(aa_aln, dna) {
  if (nchar(dna) < 3) return("")
  codons <- str_sub(dna, seq(1, nchar(dna), 3), seq(3, nchar(dna) + 2, 3))
  out <- character(0); k <- 1
  for (ch in strsplit(aa_aln, "")[[1]]) {
    if (ch == "-") out <- c(out, "---")
    else { out <- c(out, if (k <= length(codons)) codons[k] else "---"); k <- k + 1 }
  }
  paste0(out, collapse = "")
}

#' Codon-aware alignment of a DNAStringSet; returns a codon-aligned DNAStringSet.
#'
#' @param dna       a DNAStringSet (names are Genus_species or FASTA headers)
#' @param gcode     a Biostrings genetic-code vector (GENETIC_CODE or getGeneticCode("5"))
#' @param mafft_bin path to the MAFFT binary (defaults to MAFFT_BIN in the caller's env)
#' @param threads   MAFFT thread count (-1 = all cores)
codon_align <- function(dna, gcode, mafft_bin = MAFFT_BIN, threads = -1) {
  strs <- as.character(dna); names(strs) <- names(dna)
  bf  <- lapply(strs, .best_frame, gcode = gcode)
  aa  <- vapply(bf, `[[`, character(1), "aa")
  dnt <- vapply(bf, `[[`, character(1), "dna")
  keep <- nchar(aa) > 0; aa <- aa[keep]; dnt <- dnt[keep]      # drop empty translations
  aa_fa <- tempfile(fileext = ".fasta")
  writeXStringSet(AAStringSet(setNames(aa, names(aa))), aa_fa)
  aln_txt <- system2(mafft_bin, c("--auto", "--reorder", "--thread", threads, shQuote(aa_fa)),
                     stdout = TRUE, stderr = FALSE)
  aa_aln_fa <- tempfile(fileext = ".fasta"); writeLines(aln_txt, aa_aln_fa)
  aa_aln <- as.character(readAAStringSet(aa_aln_fa))
  common <- intersect(names(aa_aln), names(dnt))
  codon  <- mapply(.backtranslate, aa_aln[common], dnt[common], USE.NAMES = TRUE)
  DNAStringSet(codon)
}
