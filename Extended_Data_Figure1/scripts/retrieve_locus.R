# =====================================================
# DANAINI: RETRIEVE ONE LOCUS FROM NCBI
# - Input: species tibble with column `butterfly_species`
# - Input: one-row `locus_row` with columns: locus, term_fragment, min_len, max_len, keep_n
# - Output (per locus folder): *_accessions_raw.csv, *_accessions_trimmed.csv, *_all_trimmed.fasta
# =====================================================

options(timeout = 600)

suppressPackageStartupMessages({
  library(rentrez)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(stringr)
  library(readr)
})

`%||%` = function(x, y) if (!is.null(x) && length(x) && !is.na(x)) x else y
ncbi_pause = function() Sys.sleep(0.34)

esummary_to_list = function(sums){
  if (is.null(sums)) return(NULL)
  
  if ("extract_from_esummary" %in% getNamespaceExports("rentrez")) {
    out = tryCatch(rentrez::extract_from_esummary(sums), error = function(e) NULL)
    if (!is.null(out)) return(out)
  }
  
  if (inherits(sums, "esummary")) return(list(sums))
  if (is.list(sums)) return(sums)
  NULL
}

pick_locus_ids = function(sp,
                          term_fragment,
                          retmax = 200,
                          min_len = 200,
                          max_len = 30000){
  
  term = paste0('"', sp, '"[Organism] AND (', term_fragment, ')')
  
  s = tryCatch(
    rentrez::entrez_search(db = "nucleotide", term = term, retmax = retmax),
    error = function(e) NULL
  )
  ncbi_pause()
  
  if (is.null(s) || length(s$ids) == 0) {
    return(tibble(accession = character(), slen = integer(), n_found = 0L))
  }
  
  sums_raw = tryCatch(
    rentrez::entrez_summary(db = "nucleotide", id = s$ids),
    error = function(e) NULL
  )
  ncbi_pause()
  
  sums = esummary_to_list(sums_raw)
  if (is.null(sums)) {
    return(tibble(accession = character(), slen = integer(), n_found = as.integer(s$count)))
  }
  
  meta =
    purrr::map_dfr(sums, function(x){
      if (!is.list(x)) return(NULL)
      
      acc = x[["accessionversion"]] %||% x[["accession"]] %||% x[["uid"]]
      sl  = suppressWarnings(as.integer(x[["slen"]]))
      
      if (is.null(acc) || is.na(sl)) return(NULL)
      tibble(accession = as.character(acc), slen = sl)
    }) %>%
    filter(!is.na(accession),
           !is.na(slen),
           slen >= min_len,
           slen <= max_len) %>%
    distinct(accession, .keep_all = TRUE) %>%
    arrange(desc(slen)) %>%
    mutate(n_found = as.integer(s$count))
  
  if (nrow(meta) == 0) {
    meta = tibble(accession = character(), slen = integer(), n_found = as.integer(s$count))
  }
  
  meta
}

fetch_fasta_batch = function(ids, batch_size = 200){
  ids = unique(ids)
  if (length(ids) == 0) return(character())
  
  batches = split(ids, ceiling(seq_along(ids) / batch_size))
  
  fas = purrr::map_chr(batches, function(b){
    fa = tryCatch(
      rentrez::entrez_fetch(db = "nucleotide", id = b, rettype = "fasta", retmode = "text"),
      error = function(e) ""
    )
    ncbi_pause()
    fa
  })
  
  fas[fas != ""]
}

run_locus_pipeline = function(species,
                              locus_row,
                              retmax = 200,
                              out_dir = "ncbi_loci"){
  
  if (!("butterfly_species" %in% names(species))) {
    stop("`species` must have a column named `butterfly_species`.")
  }
  if (nrow(locus_row) != 1) {
    stop("`locus_row` must be a 1-row tibble/data.frame (filter your `loci` table to exactly one locus).")
  }
  needed = c("locus", "term_fragment", "min_len", "max_len", "keep_n")
  if (!all(needed %in% names(locus_row))) {
    stop("`locus_row` must contain columns: ", paste(needed, collapse = ", "))
  }
  
  locus = as.character(locus_row$locus[[1]])
  term_fragment = as.character(locus_row$term_fragment[[1]])
  min_len = as.integer(locus_row$min_len[[1]])
  max_len = as.integer(locus_row$max_len[[1]])
  keep_n = as.integer(locus_row$keep_n[[1]])
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  locus_dir = file.path(out_dir, locus)
  dir.create(locus_dir, showWarnings = FALSE, recursive = TRUE)
  
  message("=== ", locus, " ===")
  
  picked =
    species %>%
    mutate(
      meta = purrr::map(
        butterfly_species,
        ~ tryCatch(
          pick_locus_ids(.x,
                         term_fragment = term_fragment,
                         retmax = retmax,
                         min_len = min_len,
                         max_len = max_len),
          error = function(e) tibble(accession = character(), slen = integer(), n_found = NA_integer_)
        )
      )
    ) %>%
    mutate(
      locus_count = purrr::map_int(meta, ~ .x$n_found[1] %||% 0L),
      locus_ids   = purrr::map(meta, "accession"),
      n_ids       = purrr::map_int(locus_ids, length)
    ) %>%
    select(-meta)
  
  picked_trimmed =
    picked %>%
    mutate(
      locus_ids = purrr::map(locus_ids, head, keep_n),
      n_ids     = purrr::map_int(locus_ids, length)
    )
  
  readr::write_csv(
    picked %>%
      transmute(
        butterfly_species,
        locus_count,
        n_ids,
        accessions = purrr::map_chr(locus_ids, paste, collapse = ";")
      ),
    file.path(locus_dir, paste0(locus, "_accessions_raw.csv"))
  )
  
  readr::write_csv(
    picked_trimmed %>%
      transmute(
        butterfly_species,
        n_ids,
        accessions = purrr::map_chr(locus_ids, paste, collapse = ";")
      ),
    file.path(locus_dir, paste0(locus, "_accessions_trimmed.csv"))
  )
  
  all_acc =
    picked_trimmed %>%
    pull(locus_ids) %>%
    unlist(use.names = FALSE) %>%
    unique()
  
  message("Accessions to fetch: ", length(all_acc))
  
  fasta_all = fetch_fasta_batch(all_acc, batch_size = 200)
  fasta_path = file.path(locus_dir, paste0(locus, "_all_trimmed.fasta"))
  writeLines(fasta_all, fasta_path)
  
  miss = sum(picked_trimmed$n_ids == 0, na.rm = TRUE)
  message("Missing taxa (no hits kept): ", miss, " / ", nrow(species))
  
  invisible(list(
    locus = locus,
    locus_dir = locus_dir,
    picked = picked,
    picked_trimmed = picked_trimmed,
    fasta_path = fasta_path
  ))
}