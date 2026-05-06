# =====================================================
# drop_bad_accessions.R
# Remove listed accessions from each locus FASTA:
#   ncbi_loci/<locus>/<locus>_all_trimmed_plus_outgroups.fasta
# Writes:
#   ncbi_loci/<locus>/<locus>_all_trimmed_plus_outgroups_filtered.fasta
# Also writes a removal log:
#   curation/bad_accessions_removed_log.csv
# =====================================================

suppressPackageStartupMessages({
  library(Biostrings)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(readr)
  library(stringr)
})

# Extract accession-like tokens anywhere in a FASTA header.
# Handles: XM_032672544.2, EU141321.1, KT286348.1, OK737148.1, etc.
extract_accession = function(x){
  str_extract(x, "\\b[A-Z]{1,3}\\d{1,9}(?:_\\d+)?\\.\\d+\\b")
}

drop_bad_accessions_from_fasta = function(locus, accessions,
                                          in_dir = "ncbi_loci",
                                          in_suffix = "_all_trimmed_plus_outgroups.fasta",
                                          out_suffix = "_all_trimmed_plus_outgroups_filtered.fasta") {
  
  in_fa  = file.path(in_dir, locus, paste0(locus, in_suffix))
  out_fa = file.path(in_dir, locus, paste0(locus, out_suffix))
  
  if (!file.exists(in_fa)) {
    warning("Missing input FASTA for ", locus, ": ", in_fa)
    return(tibble(locus = locus, accession = accessions, removed = FALSE, reason = "missing_input_fasta"))
  }
  
  dna = readDNAStringSet(in_fa)
  hdr = names(dna)
  
  # Find accession anywhere in header (not just first token)
  acc = extract_accession(hdr)
  
  # If extraction fails for some headers, still allow direct string matching
  # (e.g., if headers are literally the accession as the first token)
  acc2 = ifelse(is.na(acc), sub("^([^\\s]+).*", "\\1", hdr), acc)
  
  if (!length(accessions)) {
    writeXStringSet(dna, out_fa)
    return(tibble(locus = locus, accession = character(0), removed = logical(0), reason = character(0)))
  }
  
  # Drop sequences if the extracted accession matches OR if the header contains the accession
  to_drop = acc2 %in% accessions | vapply(hdr, function(h){
    any(str_detect(h, fixed(accessions)))
  }, logical(1))
  
  dna2 = dna[!to_drop]
  writeXStringSet(dna2, out_fa)
  
  # Log which requested accessions were found (by accession token OR header contains)
  found = vapply(accessions, function(a){
    any(acc2 == a) || any(str_detect(hdr, fixed(a)))
  }, logical(1))
  
  tibble(
    locus = locus,
    accession = accessions,
    removed = found,
    reason = ifelse(found, "dropped", "not_found_in_fasta")
  )
}

run_drop_bad_accessions = function(bad_accessions,
                                   out_log = file.path("curation", "bad_accessions_removed_log.csv")) {
  
  dir.create(dirname(out_log), showWarnings = FALSE, recursive = TRUE)
  
  loci = names(bad_accessions)
  
  log_df = map_dfr(loci, ~ drop_bad_accessions_from_fasta(.x, bad_accessions[[.x]]))
  
  write_csv(log_df, out_log)
  
  invisible(log_df)
}
