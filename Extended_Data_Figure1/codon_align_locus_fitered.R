# =====================================================
# codon_align_locus_filtered.R (curation-aware)
#
# Same logic as before, but:
# - reads the FILTERED FASTA by default
# - allows choosing input suffix
# - writes outputs with a tag so you don't overwrite earlier runs
# =====================================================

suppressPackageStartupMessages({
  library(Biostrings)
  library(stringr)
})

codon_align_locus_filtered = function(locus,
                             genetic_code = Biostrings::GENETIC_CODE,
                             threads = -1,
                             in_suffix = "_all_trimmed_plus_outgroups_filtered.fasta",
                             out_tag = "filtered"){
  
  message("=== Codon-aligning ", locus, " (", out_tag, ") ===")
  
  # -------- paths --------
  in_fa = file.path("ncbi_loci", locus, paste0(locus, in_suffix))
  if (!file.exists(in_fa)) stop("Input FASTA not found: ", in_fa)
  
  outdir = file.path("ncbi_alignments_filtered", locus)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  aa_fa     = file.path(outdir, paste0(locus, "_", out_tag, "_aa.fasta"))
  aa_aln_fa = file.path(outdir, paste0(locus, "_", out_tag, "_aa_aln.fasta"))
  codon_aln = file.path(outdir, paste0(locus, "_", out_tag, "_codon_aln.fasta"))
  
  trim_to_full_codons = function(s){
    s = toupper(gsub("\\s+", "", s))
    n3 = nchar(s) - (nchar(s) %% 3)
    if (n3 < 3) return("")
    substr(s, 1, n3)
  }
  
  translate_frame = function(dna_str, frame = 1){
    dna_str = toupper(gsub("\\s+", "", dna_str))
    if (nchar(dna_str) < frame) return("")
    s = substr(dna_str, frame, nchar(dna_str))
    s = trim_to_full_codons(s)
    if (nchar(s) < 3) return("")
    aa = suppressWarnings(
      Biostrings::translate(
        DNAString(s),
        genetic.code = genetic_code,
        if.fuzzy.codon = "X"
      )
    )
    as.character(aa)
  }
  
  score_translation = function(aa){
    if (nchar(aa) == 0) return(list(stops = Inf, len = 0))
    stops = stringr::str_count(aa, fixed("*"))
    if (endsWith(aa, "*")) stops = stops - 1
    list(stops = stops, len = nchar(aa))
  }
  
  best_frame = function(dna_str){
    a1 = translate_frame(dna_str, 1)
    a2 = translate_frame(dna_str, 2)
    a3 = translate_frame(dna_str, 3)
    
    s1 = score_translation(a1)
    s2 = score_translation(a2)
    s3 = score_translation(a3)
    
    df = data.frame(
      frame = c(1,2,3),
      stops = c(s1$stops, s2$stops, s3$stops),
      len   = c(s1$len,   s2$len,   s3$len)
    )
    df = df[order(df$stops, -df$len), ]
    best = df$frame[1]
    
    dna_in_frame = substr(toupper(dna_str), best, nchar(dna_str))
    dna_in_frame = trim_to_full_codons(dna_in_frame)
    aa_best = list(a1,a2,a3)[[best]]
    
    list(aa = aa_best, dna = dna_in_frame)
  }
  
  backtranslate = function(aa_aln, dna){
    if (nchar(dna) < 3) return("")
    starts = seq(1, nchar(dna), by = 3)
    codons = stringr::str_sub(dna, starts, starts + 2)
    
    out = character(0)
    k = 1
    for (ch in strsplit(aa_aln, "")[[1]]){
      if (ch == "-"){
        out = c(out, "---")
      } else {
        if (k <= length(codons)){
          out = c(out, codons[k]); k = k + 1
        } else {
          out = c(out, "---")
        }
      }
    }
    paste0(out, collapse = "")
  }
  
  # 1) read nucleotide sequences
  dna = readDNAStringSet(in_fa)
  dna_strs = as.character(dna)
  names(dna_strs) = names(dna)
  
  # 2) choose best frame per sequence + translate
  bf = lapply(dna_strs, best_frame)
  aa_strs  = vapply(bf, `[[`, character(1), "aa")
  dna_trim = vapply(bf, `[[`, character(1), "dna")
  
  # drop empty translations (important for MAFFT)
  keep = nchar(aa_strs) > 0
  aa_strs  = aa_strs[keep]
  dna_trim = dna_trim[keep]
  
  aa_set = AAStringSet(aa_strs)
  names(aa_set) = names(aa_strs)
  writeXStringSet(aa_set, aa_fa)
  
  # 3) align amino acids with MAFFT
  if (!nzchar(Sys.which("mafft"))) stop("MAFFT not found in PATH.")
  message("Running MAFFT...")
  
  aln_txt = system2(
    "mafft",
    c("--auto", "--reorder", "--thread", threads, shQuote(aa_fa)),
    stdout = TRUE
  )
  writeLines(aln_txt, aa_aln_fa)
  
  aa_aln = readAAStringSet(aa_aln_fa)
  aa_aln_str = as.character(aa_aln)
  
  # 4) back-translate
  common = intersect(names(aa_aln_str), names(dna_trim))
  
  codon_str = mapply(
    backtranslate,
    aa_aln_str[common],
    dna_trim[common],
    USE.NAMES = TRUE
  )
  
  codon_set = DNAStringSet(codon_str)
  names(codon_set) = names(codon_str)
  writeXStringSet(codon_set, codon_aln)
  
  message("✓ Finished ", locus)
  invisible(list(aa_aln = aa_aln_fa, codon_aln = codon_aln))
}