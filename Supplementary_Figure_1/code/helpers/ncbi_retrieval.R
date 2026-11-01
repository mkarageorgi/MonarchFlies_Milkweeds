#!/usr/bin/env Rscript
# =====================================================
# ncbi_retrieval.R
#
# Retrieve sequence data for a locus from NCBI nucleotide, for a set of butterfly species.
#
#   run_locus_pipeline(species, locus_row, out_dir)   generic nuclear locus
#   run_coi_pipeline(species, COI_row, out_dir)        COI (mitogenome-safe + COII-safe)
#   get_danaini_species()                              Danaini species list from NCBI taxonomy
#
# Per-locus outputs (in out_dir/<locus>/):
#   <locus>_accessions_raw.csv, <locus>_accessions_trimmed.csv, <locus>_all_trimmed.fasta
#
# After every batch of NCBI calls, inspect the error log written by the caller: a species that
# FAILED a request (n_found = NA, logged) is not the same as one with no sequence (n_found = 0).
#
# Provenance:
#   The retrieval core is Marianna Karageorgi's, originally split across retrieve_locus.R
#   (run_locus_pipeline, pick_locus_ids, fetch_fasta_batch, esummary_to_list, %||%) and
#   retrieve_coi.R (run_coi_pipeline, pick_coi_ids). N. Weksler added the retry + error-logging
#   layer (ncbi_try / log_ncbi_error, and the robust pick_*/fetch_* that wrap every bare entrez_*
#   call) and the species-sampling functions (get_danaini_taxid / get_danaini_species), so that a
#   transient network/throttle failure is distinguishable from a genuine absence.
#
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

`%||%` <- function(x, y) if (!is.null(x) && length(x) && !is.na(x)) x else y

# Faster pacing is safe with an ENTREZ_KEY registered (~10 req/s).
ncbi_pause <- function() Sys.sleep(0.11)

esummary_to_list <- function(sums) {
  if (is.null(sums)) return(NULL)
  if ("extract_from_esummary" %in% getNamespaceExports("rentrez")) {
    out <- tryCatch(rentrez::extract_from_esummary(sums), error = function(e) NULL)
    if (!is.null(out)) return(out)
  }
  if (inherits(sums, "esummary")) return(list(sums))
  if (is.list(sums)) return(sums)
  NULL
}

# ---- retry + error log --------------------------------------------------------
# One row per NCBI call that failed after all retries.
.ncbi_err <- new.env(parent = emptyenv())
.ncbi_err$rows <- list()
log_ncbi_error <- function(context, msg) {
  .ncbi_err$rows[[length(.ncbi_err$rows) + 1L]] <- tibble(
    context = context, message = substr(msg, 1, 300), time = as.character(Sys.time())
  )
}
# Collect the accumulated error log as a tibble (empty tibble if no failures).
ncbi_error_log <- function() {
  if (length(.ncbi_err$rows)) bind_rows(.ncbi_err$rows) else
    tibble(context = character(), message = character(), time = character())
}

#' Run one NCBI call with retries + exponential backoff.
#' @param fun zero-arg function performing a single entrez_* call
#' @return the call's result, or NULL if it errored on every attempt (then logged)
ncbi_try <- function(fun, context = "", max_tries = 3L, base_wait = 0.6) {
  for (attempt in seq_len(max_tries)) {
    res <- tryCatch(fun(), error = function(e) structure(conditionMessage(e), class = "ncbi_fail"))
    if (!inherits(res, "ncbi_fail")) return(res)
    if (attempt < max_tries) Sys.sleep(base_wait * 2^(attempt - 1L))
    else { log_ncbi_error(context, as.character(res)); return(NULL) }
  }
}

# ---- robust pickers / fetcher (every entrez_* call wrapped in ncbi_try) --------
# Sentinel: n_found = NA marks a retrieval ERROR for that species x locus (vs 0 = genuine absence).

pick_locus_ids <- function(sp, term_fragment, retmax = 200, min_len = 200, max_len = 30000) {
  term <- paste0('"', sp, '"[Organism] AND (', term_fragment, ')')
  s <- ncbi_try(function() rentrez::entrez_search(db = "nucleotide", term = term, retmax = retmax),
                context = paste0("search|", sp, "|", substr(term_fragment, 1, 24)))
  ncbi_pause()
  if (is.null(s))         return(tibble(accession = character(), slen = integer(), n_found = NA_integer_))
  if (length(s$ids) == 0) return(tibble(accession = character(), slen = integer(), n_found = 0L))

  sums_raw <- ncbi_try(function() rentrez::entrez_summary(db = "nucleotide", id = s$ids),
                       context = paste0("summary|", sp))
  ncbi_pause()
  sums <- esummary_to_list(sums_raw)
  if (is.null(sums)) return(tibble(accession = character(), slen = integer(), n_found = as.integer(s$count)))

  meta <- purrr::map_dfr(sums, function(x) {
    if (!is.list(x)) return(NULL)
    acc <- x[["accessionversion"]] %||% x[["accession"]] %||% x[["uid"]]
    sl  <- suppressWarnings(as.integer(x[["slen"]]))
    if (is.null(acc) || is.na(sl)) return(NULL)
    tibble(accession = as.character(acc), slen = sl)
  }) %>%
    filter(!is.na(accession), !is.na(slen), slen >= min_len, slen <= max_len) %>%
    distinct(accession, .keep_all = TRUE) %>%
    arrange(desc(slen)) %>%
    mutate(n_found = as.integer(s$count))

  if (nrow(meta) == 0) meta <- tibble(accession = character(), slen = integer(), n_found = as.integer(s$count))
  meta
}

pick_coi_ids <- function(sp, term_fragment, retmax = 200, min_len = 200, max_len = 30000) {
  sp_safe <- gsub('"', '\\"', sp)
  term <- paste0(
    '"', sp_safe, '"[Organism] AND (', term_fragment, ') ',
    'NOT ("mitochondrion, complete genome"[Title] OR "mitochondrial genome"[Title] OR ',
    '"complete mitochondrion"[Title] OR "complete genome"[Title] OR WGS[Keyword]) ',
    'NOT (COII[Title] OR COX2[Title] OR cox2[Title] OR "cytochrome c oxidase subunit II"[Title] OR ',
    '"cytochrome oxidase subunit II"[Title])'
  )
  s <- ncbi_try(function() rentrez::entrez_search(db = "nucleotide", term = term, retmax = retmax),
                context = paste0("search-COI|", sp))
  ncbi_pause()
  if (is.null(s))         return(tibble(accession = character(), slen = integer(), n_found = NA_integer_))
  if (length(s$ids) == 0) return(tibble(accession = character(), slen = integer(), n_found = 0L))

  sums_raw <- ncbi_try(function() rentrez::entrez_summary(db = "nucleotide", id = s$ids),
                       context = paste0("summary-COI|", sp))
  ncbi_pause()
  sums <- esummary_to_list(sums_raw)
  if (is.null(sums)) return(tibble(accession = character(), slen = integer(), n_found = as.integer(s$count)))

  meta <- purrr::map_dfr(sums, function(x) {
    if (!is.list(x)) return(NULL)
    acc <- x[["accessionversion"]] %||% x[["accession"]] %||% x[["uid"]]
    sl  <- suppressWarnings(as.integer(x[["slen"]]))
    ttl <- x[["title"]] %||% ""
    if (is.null(acc) || is.na(sl)) return(NULL)
    tibble(accession = as.character(acc), slen = sl, title = as.character(ttl))
  }) %>%
    filter(!is.na(accession), !is.na(slen), slen >= min_len, slen <= max_len) %>%
    filter(!str_detect(str_to_lower(title),
      "mitochondrion, complete genome|mitochondrial genome|complete mitochondrion|genome, complete|complete genome|whole genome")) %>%
    filter(!str_detect(str_to_lower(title),
      "\\bcoii\\b|\\bcox2\\b|cytochrome c oxidase subunit ii|cytochrome oxidase subunit ii")) %>%
    distinct(accession, .keep_all = TRUE) %>%
    arrange(desc(slen)) %>%
    mutate(n_found = as.integer(s$count)) %>%
    select(accession, slen, n_found)

  if (nrow(meta) == 0) meta <- tibble(accession = character(), slen = integer(), n_found = as.integer(s$count))
  meta
}

fetch_fasta_batch <- function(ids, batch_size = 200) {
  ids <- unique(ids)
  if (length(ids) == 0) return(character())
  batches <- split(ids, ceiling(seq_along(ids) / batch_size))
  fas <- purrr::map_chr(batches, function(b) {
    fa <- ncbi_try(function() rentrez::entrez_fetch(db = "nucleotide", id = b,
                                                    rettype = "fasta", retmode = "text"),
                   context = paste0("fetch|", length(b), " ids"))
    ncbi_pause()
    if (is.null(fa)) "" else fa
  })
  fas[fas != ""]
}

# ---- per-locus pipelines (call the pickers/fetcher above by bare name) ---------

run_locus_pipeline <- function(species, locus_row, retmax = 200, out_dir = "ncbi_loci") {
  if (!("butterfly_species" %in% names(species)))
    stop("`species` must have a column named `butterfly_species`.")
  if (nrow(locus_row) != 1)
    stop("`locus_row` must be a 1-row tibble (filter your `loci` table to exactly one locus).")
  needed <- c("locus", "term_fragment", "min_len", "max_len", "keep_n")
  if (!all(needed %in% names(locus_row)))
    stop("`locus_row` must contain columns: ", paste(needed, collapse = ", "))

  locus         <- as.character(locus_row$locus[[1]])
  term_fragment <- as.character(locus_row$term_fragment[[1]])
  min_len <- as.integer(locus_row$min_len[[1]]); max_len <- as.integer(locus_row$max_len[[1]])
  keep_n  <- as.integer(locus_row$keep_n[[1]])

  locus_dir <- file.path(out_dir, locus); dir.create(locus_dir, showWarnings = FALSE, recursive = TRUE)
  message("=== ", locus, " ===")

  picked <- species %>%
    mutate(meta = purrr::map(butterfly_species, ~ tryCatch(
      pick_locus_ids(.x, term_fragment, retmax, min_len, max_len),
      error = function(e) tibble(accession = character(), slen = integer(), n_found = NA_integer_)))) %>%
    mutate(locus_count = purrr::map_int(meta, ~ .x$n_found[1] %||% 0L),
           locus_ids   = purrr::map(meta, "accession"),
           n_ids       = purrr::map_int(locus_ids, length)) %>%
    select(-meta)

  picked_trimmed <- picked %>%
    mutate(locus_ids = purrr::map(locus_ids, head, keep_n),
           n_ids     = purrr::map_int(locus_ids, length))

  write_csv(picked %>% transmute(butterfly_species, locus_count, n_ids,
              accessions = purrr::map_chr(locus_ids, paste, collapse = ";")),
            file.path(locus_dir, paste0(locus, "_accessions_raw.csv")))
  write_csv(picked_trimmed %>% transmute(butterfly_species, n_ids,
              accessions = purrr::map_chr(locus_ids, paste, collapse = ";")),
            file.path(locus_dir, paste0(locus, "_accessions_trimmed.csv")))

  all_acc <- picked_trimmed %>% pull(locus_ids) %>% unlist(use.names = FALSE) %>% unique()
  message("Accessions to fetch: ", length(all_acc))
  fasta_path <- file.path(locus_dir, paste0(locus, "_all_trimmed.fasta"))
  writeLines(fetch_fasta_batch(all_acc, batch_size = 200), fasta_path)
  message("Missing species (no hits kept): ", sum(picked_trimmed$n_ids == 0, na.rm = TRUE),
          " / ", nrow(species))

  invisible(list(locus = locus, locus_dir = locus_dir, picked = picked,
                 picked_trimmed = picked_trimmed, fasta_path = fasta_path))
}

run_coi_pipeline <- function(species, COI_row, retmax = 200, out_dir = "ncbi_loci") {
  if (!("butterfly_species" %in% names(species)))
    stop("`species` must have a column named `butterfly_species`.")
  if (nrow(COI_row) != 1)
    stop("`COI_row` must be a 1-row tibble (filter your `loci` table to exactly one locus).")
  needed <- c("locus", "term_fragment", "min_len", "max_len", "keep_n")
  if (!all(needed %in% names(COI_row)))
    stop("`COI_row` must contain columns: ", paste(needed, collapse = ", "))

  locus <- as.character(COI_row$locus[[1]])
  if (locus != "COI") stop("run_coi_pipeline() expects COI_row$locus == 'COI'.")
  term_fragment <- as.character(COI_row$term_fragment[[1]])
  min_len <- max(as.integer(COI_row$min_len[[1]]), 500L)   # COI: enforce safe length bounds
  max_len <- min(as.integer(COI_row$max_len[[1]]), 2500L)
  keep_n  <- as.integer(COI_row$keep_n[[1]])

  locus_dir <- file.path(out_dir, locus); dir.create(locus_dir, showWarnings = FALSE, recursive = TRUE)
  message("=== ", locus, " (COI-safe, no COII) ===")

  picked <- species %>%
    mutate(meta = purrr::map(butterfly_species, ~ tryCatch(
      pick_coi_ids(.x, term_fragment, retmax, min_len, max_len),
      error = function(e) tibble(accession = character(), slen = integer(), n_found = NA_integer_)))) %>%
    mutate(locus_count = purrr::map_int(meta, ~ .x$n_found[1] %||% 0L),
           locus_ids   = purrr::map(meta, "accession"),
           n_ids       = purrr::map_int(locus_ids, length)) %>%
    select(-meta)

  picked_trimmed <- picked %>%
    mutate(locus_ids = purrr::map(locus_ids, head, keep_n),
           n_ids     = purrr::map_int(locus_ids, length))

  write_csv(picked %>% transmute(butterfly_species, locus_count, n_ids,
              accessions = purrr::map_chr(locus_ids, paste, collapse = ";")),
            file.path(locus_dir, paste0(locus, "_accessions_raw.csv")))
  write_csv(picked_trimmed %>% transmute(butterfly_species, n_ids,
              accessions = purrr::map_chr(locus_ids, paste, collapse = ";")),
            file.path(locus_dir, paste0(locus, "_accessions_trimmed.csv")))

  all_acc <- picked_trimmed %>% pull(locus_ids) %>% unlist(use.names = FALSE) %>% unique()
  message("Accessions to fetch: ", length(all_acc))
  fasta_path <- file.path(locus_dir, paste0(locus, "_all_trimmed.fasta"))
  writeLines(fetch_fasta_batch(all_acc, batch_size = 200), fasta_path)
  message("Missing species (no hits kept): ", sum(picked_trimmed$n_ids == 0, na.rm = TRUE),
          " / ", nrow(species))

  invisible(list(locus = locus, locus_dir = locus_dir, picked = picked,
                 picked_trimmed = picked_trimmed, fasta_path = fasta_path))
}

# ---- species sampling from NCBI taxonomy --------------------------------------

get_danaini_taxid <- function() {
  hit <- ncbi_try(function() rentrez::entrez_search(db = "taxonomy",
                    term = "Danaini[Scientific Name]", retmax = 10), context = "taxid-search")
  ncbi_pause()
  if (is.null(hit)) stop("Could not resolve Danaini taxid (NCBI error).")
  sums <- ncbi_try(function() rentrez::entrez_summary(db = "taxonomy", id = hit$ids),
                   context = "taxid-summary")
  ncbi_pause()
  sums <- if (inherits(sums, "esummary")) list(sums) else sums
  tribe <- purrr::keep(sums, ~ tolower(.x$rank %||% "") == "tribe" &&
                              (.x$scientificname %||% "") == "Danaini")
  as.character((tribe[[1]]$uid %||% hit$ids[[1]]))
}

get_danaini_species <- function(retmax = 20000) {
  danaini_id <- get_danaini_taxid()
  desc <- ncbi_try(function() rentrez::entrez_search(db = "taxonomy",
                     term = paste0("txid", danaini_id, "[Subtree]"), retmax = retmax),
                   context = "danaini-subtree")
  ncbi_pause()
  if (is.null(desc) || length(desc$ids) == 0) stop("No descendants returned for Danaini subtree.")

  batches <- split(desc$ids, ceiling(seq_along(desc$ids) / 500))
  sums <- purrr::map(batches, function(b) {
    out <- ncbi_try(function() rentrez::entrez_summary(db = "taxonomy", id = b),
                    context = "danaini-tax-summary")
    ncbi_pause()
    if (inherits(out, "esummary")) out <- list(out)
    out
  }) %>% unlist(recursive = FALSE)

  tibble(
    butterfly_species = map_chr(sums, ~ .x$scientificname %||% NA_character_),
    rank              = map_chr(sums, ~ tolower(.x$rank %||% ""))
  ) %>%
    filter(rank == "species") %>%
    transmute(butterfly_species) %>%
    filter(str_detect(butterfly_species, "^[A-Z][a-z]+\\s+[a-z]+$")) %>%
    distinct(butterfly_species) %>%
    arrange(butterfly_species)
}
