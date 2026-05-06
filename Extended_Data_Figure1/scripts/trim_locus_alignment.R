# =====================================================
# trim_codon_alignment.R
# Trim codon-aligned FASTA to a start/end window, filter gappy seqs,
# optionally collapse to one seq per species, and write clean headers:
#   >Genus_species
#
# Also writes an accession lookup table alongside each output FASTA:
#   <out_fasta>_accession_map.csv
# =====================================================

suppressPackageStartupMessages({
  library(Biostrings)
  library(stringr)
})

# ---- parse accession token from header (best-effort) ----
get_accession = function(hdr){
  # tries to extract GenBank-like accession anywhere; falls back to first token
  acc = str_extract(hdr, "\\b[A-Z]{1,3}\\d{1,9}(?:_\\d+)?\\.\\d+\\b")
  ifelse(is.na(acc), sub("^([^\\s]+).*", "\\1", hdr), acc)
}

# ---- extract binomial from header and convert to Genus_species ----
standardize_taxon = function(hdr){
  sp = str_extract(hdr, "[A-Z][a-z]+\\s+[a-z]+")
  if (is.na(sp)) return(NA_character_)
  gsub("\\s+", "_", sp)
}

# ---- fraction of gap/ambig ----
bad_fraction = function(x){
  s = strsplit(as.character(x), "", fixed = TRUE)[[1]]
  mean(s %in% c("-", "N", "n", "?"))
}

trim_codon_alignment = function(in_fasta,
                                out_fasta,
                                start,
                                end,
                                max_bad = 0.20,
                                one_per_species = TRUE,
                                standardize_names = TRUE,
                                fail_on_unknown = TRUE,
                                write_accession_map = TRUE){
  
  # 0) codon-safe trimming
  if (((start - 1) %% 3) != 0)
    stop("start must satisfy (start - 1) %% 3 == 0 for codon alignment (got start=", start, ")")
  if (((end - start + 1) %% 3) != 0)
    stop("(end - start + 1) must be divisible by 3 (got start=", start, ", end=", end, ")")
  if (end < start) stop("end < start")
  
  # 1) read
  aln = readDNAStringSet(in_fasta)
  
  # 2) trim columns
  aln_trim = subseq(aln, start = start, end = end)
  
  # 3) gap/ambig filter
  frac_bad = vapply(aln_trim, bad_fraction, numeric(1))
  keep_gaps = frac_bad <= max_bad
  aln_filt = aln_trim[keep_gaps]
  frac_bad = frac_bad[keep_gaps]
  
  # 4) standardize names (for cross-locus concatenation)
  hdr2 = names(aln_filt)
  acc2 = get_accession(hdr2)
  tax  = vapply(hdr2, standardize_taxon, character(1))
  
  if (standardize_names) {
    bad = which(is.na(tax) | tax == "")
    if (length(bad) > 0) {
      msg = paste0("Could not parse binomial from ", length(bad),
                   " headers (first 5 shown):\n  - ",
                   paste(head(hdr2[bad], 5), collapse = "\n  - "))
      if (fail_on_unknown) stop(msg) else warning(msg)
    }
  }
  
  # 5) one sequence per species (choose best by lowest gap, then longest)
  if (one_per_species) {
    
    sp_key = ifelse(is.na(tax), "UNKNOWN", tax)
    seq_len = width(aln_filt)
    gap = frac_bad
    
    df = data.frame(
      idx = seq_along(aln_filt),
      species = sp_key,
      accession = acc2,
      gap = gap,
      len = seq_len,
      stringsAsFactors = FALSE
    )
    
    df = df[order(df$species, df$gap, -df$len), ]
    best = !duplicated(df$species)
    
    aln_filt = aln_filt[df$idx[best]]
    df_best = df[best, , drop = FALSE]
    
    message("After species collapsing: ", length(aln_filt), " sequences (one per species)")
    
  } else {
    df_best = data.frame(
      species = ifelse(is.na(tax), "UNKNOWN", tax),
      accession = acc2,
      gap = frac_bad,
      len = width(aln_filt),
      stringsAsFactors = FALSE
    )
  }
  
  # 6) rename headers cleanly: >Genus_species
  if (standardize_names) {
    clean_names = ifelse(is.na(df_best$species) | df_best$species == "", "UNKNOWN", df_best$species)
    names(aln_filt) = clean_names
  } else {
    names(aln_filt) = sub("^(\\S+).*$", "\\1", names(aln_filt))
  }
  
  # 7) write FASTA
  dir.create(dirname(out_fasta), showWarnings = FALSE, recursive = TRUE)
  writeXStringSet(aln_filt, out_fasta)
  
  # 8) write accession map
  if (write_accession_map) {
    map_path = sub("\\.fasta$", "_accession_map.csv", out_fasta)
    if (identical(map_path, out_fasta)) map_path = paste0(out_fasta, "_accession_map.csv")
    write.csv(df_best[, c("species","accession","gap","len")], map_path, row.names = FALSE)
  }
  
  message("Wrote: ", out_fasta)
  message("  n_in: ", length(aln),
          " | n_after_gapheavy: ", sum(keep_gaps),
          " | final (after 1/sp): ", length(aln_filt),
          " | width: ", unique(width(aln_filt)))
  
  invisible(list(
    n_in = length(aln),
    n_after_gapheavy = sum(keep_gaps),
    n_final = length(aln_filt),
    width = unique(width(aln_filt)),
    removed_gapheavy = sum(!keep_gaps)
  ))
}
