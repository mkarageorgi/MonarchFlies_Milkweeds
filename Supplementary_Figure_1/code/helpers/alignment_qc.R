#!/usr/bin/env Rscript
# =====================================================
# alignment_qc.R
#
# Automated trimming + quality control of a codon alignment, and a per-locus convenience
# wrapper that runs read -> align -> QC for one locus.
#
#   auto_trim_and_qc(aln, gcode, ...) -> list(aln = <QC'd DNAStringSet>, qc = <log>, n_cols)
#   process_locus(locus, keep_std, ...) -> QC'd DNAStringSet for one locus (or NULL)
#
# What the QC does:
#   - trims to the dense codon-column core (columns present in >= occ_threshold of sequences),
#     preserving reading frame and stripping MAFFT's sparse insertions / mRNA-only tails;
#   - drops a sequence only if it is mostly gap over that core (max_bad), has >1 internal stop,
#     or is a consensus distance-outlier;
#   - protects outgroup species from every drop filter (they are legitimately divergent and are
#     needed to root the tree);
#   - keeps one sequence per species (lowest gap fraction, then most non-gap bases).
#
# Provenance:
#   Written by N. Weksler. This automated QC replaces the manual "inspect the alignment in
#   Galaxy and list bad accessions by hand" step in Marianna Karageorgi's original workflow;
#   every dropped sequence and the reason for dropping it is now recorded in a per-locus
#   <locus>_qc_log.csv so the thresholds can be tuned objectively.
#
# Depends on helpers in sequence_utils.R (standardize_taxon, get_accession, read_combined_locus,
# subset_to_species) and codon_alignment.R (codon_align) being sourced first.
#
# =====================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(Biostrings)
})

# Codon-column occupancy trim: keep only codon columns (triplets) that are fully present
# (all 3 positions non-gap) in >= occ_threshold of sequences. Operating in codon units
# preserves reading frame and removes MAFFT's sparse insertions and long mRNA-only tails
# that would otherwise make short amplicons look gappy and get dropped.
.keep_codon_cols <- function(mat, occ_threshold) {
  ncodon <- ncol(mat) %/% 3
  keep <- vapply(seq_len(ncodon), function(c) {
    cols <- (3 * c - 2):(3 * c)
    mean(rowSums(mat[, cols, drop = FALSE] == "-") == 0) >= occ_threshold
  }, logical(1))
  if (!any(keep)) keep[] <- TRUE   # fallback: keep all codons if none pass
  as.integer(unlist(lapply(which(keep), function(c) (3 * c - 2):(3 * c))))
}

# Count internal stop codons in a codon-aligned nucleotide string. `gcode` is a named
# codon->AA vector, so this is a vectorized lookup (not a per-codon translate() call, which
# was ~1000x slower).
.count_internal_stops <- function(seq_chr, gcode) {
  cod <- substring(seq_chr, seq(1, nchar(seq_chr) - 2, 3), seq(3, nchar(seq_chr), 3))
  aa  <- gcode[toupper(cod)]          # NA for gap / ambiguous codons
  aa[is.na(aa)] <- "X"
  n <- sum(aa == "*")
  if (length(aa) && aa[length(aa)] == "*") n <- n - 1L   # ignore a terminal stop
  as.integer(n)
}

# p-distance of each sequence to the column-majority consensus; flag outliers beyond
# median + mad_k * MAD.
.consensus_outliers <- function(mat, mad_k) {
  cons <- apply(mat, 2, function(col) {
    v <- col[!(col %in% c("-", "N", "n", "?"))]
    if (!length(v)) "-" else names(sort(table(v), decreasing = TRUE))[1]
  })
  d <- apply(mat, 1, function(row) {
    ok <- !(row %in% c("-", "N", "n", "?")) & cons != "-"
    if (!any(ok)) NA_real_ else mean(row[ok] != cons[ok])
  })
  thr <- median(d, na.rm = TRUE) + mad_k * mad(d, na.rm = TRUE)
  list(dist = d, threshold = thr, outlier = !is.na(d) & d > thr)
}

#' Trim + QC a codon alignment.
#' @return list(aln = <QC'd DNAStringSet, names = Genus_species>, qc = <tibble log>, n_cols)
auto_trim_and_qc <- function(aln, gcode, protect_species = character(0),
                             max_bad = MAX_BAD, occ_threshold = OCC_THRESHOLD,
                             mad_k = MAD_K, locus = "") {
  hdr <- names(aln)
  species <- vapply(hdr, standardize_taxon, character(1))
  acc     <- vapply(hdr, get_accession, character(1))
  ok <- !is.na(species)
  species <- species[ok]; acc <- acc[ok]
  mat <- do.call(rbind, strsplit(as.character(aln[ok]), "", fixed = TRUE))  # nseq x ncol

  keep_cols <- .keep_codon_cols(mat, occ_threshold)
  mat_t <- mat[, keep_cols, drop = FALSE]
  seqstr <- apply(mat_t, 1, paste0, collapse = "")

  protected <- species %in% protect_species
  badf  <- rowMeans(matrix(mat_t %in% c("-", "N", "n", "?"), nrow = nrow(mat_t)))
  stops <- vapply(seqstr, .count_internal_stops, integer(1), gcode = gcode)
  co    <- .consensus_outliers(mat_t, mad_k)

  drop <- (badf > max_bad) | (stops > 1) | co$outlier
  drop[protected] <- FALSE

  # one sequence per species: lowest gap fraction, then most non-gap bases
  keep <- !drop
  sel_idx <- integer(0)                        # indices (into the parsed-species set) chosen
  if (any(keep)) {
    kept_orig <- which(keep)
    sp_k <- species[keep]; bad_k <- badf[keep]; seq_k <- seqstr[keep]
    len_k <- nchar(gsub("-", "", seq_k))
    ord  <- order(sp_k, bad_k, -len_k)
    sel_idx <- kept_orig[ord[!duplicated(sp_k[ord])]]
  }
  selected <- logical(length(species)); selected[sel_idx] <- TRUE

  qc_log <- tibble(locus = locus, species = species, accession = acc,
                   bad_frac = round(badf, 3), internal_stops = stops,
                   cons_dist = round(co$dist, 3), cons_threshold = round(co$threshold, 3),
                   protected = protected, dropped = drop,
                   reason = case_when(
                     !drop ~ "",
                     badf > max_bad ~ "gap/ambiguity",
                     stops > 1 ~ "internal_stops",
                     co$outlier ~ "consensus_outlier",
                     TRUE ~ "other"),
                   selected = selected)   # TRUE = the accession that entered the supermatrix

  aln_f <- DNAStringSet(seqstr[sel_idx]); names(aln_f) <- species[sel_idx]
  list(aln = aln_f, qc = qc_log, n_cols = length(keep_cols))
}

#' Full per-locus pipeline (read -> codon-align -> QC); writes the per-locus QC log.
#' @param curation_dir where to write <locus>_qc_log.csv (default CURATION_DIR in caller's env)
process_locus <- function(locus, keep_std = NULL, curation_dir = CURATION_DIR) {
  dna <- read_combined_locus(locus)
  if (!is.null(keep_std)) dna <- subset_to_species(dna, keep_std)
  if (length(dna) < 4) { message("  <4 sequences for ", locus, " — skipped"); return(NULL) }
  aln <- codon_align(dna, GCODE_OF(locus))
  res <- auto_trim_and_qc(aln, GCODE_OF(locus), protect_species = OUTGROUPS_STD, locus = locus)
  write_csv(res$qc, file.path(curation_dir, paste0(locus, "_qc_log.csv")))
  if (length(res$aln) < 4) { message("  ", locus, ": <4 species after QC — skipped"); return(NULL) }
  message("  ", locus, ": ", length(res$aln), " species | core ", res$n_cols,
          " bp | dropped ", sum(res$qc$dropped))
  res$aln
}
