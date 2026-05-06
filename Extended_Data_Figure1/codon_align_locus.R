# =====================================================
# codon_align_locus.R
#
# Codon-aware alignment for protein-coding loci:
#   1) read nucleotide sequences (FASTA)
#   2) pick best reading frame per sequence (1/2/3) by:
#        - minimizing internal stop codons
#        - maximizing translated length
#   3) translate to amino acids
#   4) align amino acids with MAFFT (more reliable than nucleotide alignment)
#   5) back-translate AA alignment to a codon-aligned nucleotide alignment
#
# Output per locus:
#   ncbi_alignments/<locus>/<locus>_aa_aln.fasta
#   ncbi_alignments/<locus>/<locus>_codon_aln.fasta
#
# Notes:
# - For nuclear loci, use genetic_code = Biostrings::GENETIC_CODE (default)
# - For COI (invertebrate mitochondrion), use genetic_code = Biostrings::getGeneticCode("5")
# - This code avoids width()/subseq() so it works even if S4 width methods misbehave.
# =====================================================

suppressPackageStartupMessages({
  library(Biostrings)  # DNA/AA sequence I/O + translate()
  library(stringr)     # string helpers (counting stops, extracting codons)
})

# -----------------------------------------------------
# Main function: run codon-aware alignment for ONE locus
# -----------------------------------------------------
codon_align_locus = function(locus,
                             genetic_code = Biostrings::GENETIC_CODE, # nuclear default
                             threads = -1){                             # MAFFT: -1 uses all cores
  
  message("=== Codon-aligning ", locus, " ===")
  
  # -------- paths --------
  # Input is your combined FASTA (Danaini + outgroups) produced earlier
  in_fa  = file.path("ncbi_loci", locus,
                     paste0(locus, "_all_trimmed_plus_outgroups.fasta"))
  
  # Output directory for this locus
  outdir = file.path("ncbi_alignments", locus)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  # Intermediate/output files
  aa_fa     = file.path(outdir, paste0(locus, "_aa.fasta"))       # translated AA sequences (unaligned)
  aa_aln_fa = file.path(outdir, paste0(locus, "_aa_aln.fasta"))   # AA alignment (MAFFT output)
  codon_aln = file.path(outdir, paste0(locus, "_codon_aln.fasta"))# codon-aligned nucleotide alignment
  
  # Basic input sanity check
  if (!file.exists(in_fa))
    stop("Input FASTA not found: ", in_fa)
  
  # -----------------------------------------------------
  # Helper: trim a DNA string to a whole number of codons
  # -----------------------------------------------------
  trim_to_full_codons = function(s){
    # remove whitespace, uppercase
    s = toupper(gsub("\\s+", "", s))
    # keep only the largest prefix whose length is divisible by 3
    n3 = nchar(s) - (nchar(s) %% 3)
    if (n3 < 3) return("")
    substr(s, 1, n3)
  }
  
  # -----------------------------------------------------
  # Helper: translate a DNA string in a given frame (1/2/3)
  # -----------------------------------------------------
  translate_frame = function(dna_str, frame = 1){
    dna_str = toupper(gsub("\\s+", "", dna_str))
    if (nchar(dna_str) < frame) return("")  # too short for this frame
    
    # shift into frame and trim to full codons
    s = substr(dna_str, frame, nchar(dna_str))
    s = trim_to_full_codons(s)
    if (nchar(s) < 3) return("")
    
    # translate; ambiguous/gappy codons become 'X'
    aa = suppressWarnings(
      Biostrings::translate(
        DNAString(s),
        genetic.code = genetic_code,
        if.fuzzy.codon = "X"
      )
    )
    as.character(aa)
  }
  
  # -----------------------------------------------------
  # Helper: score a translation
  #   - fewer internal stop codons is better
  #   - longer AA sequence is better
  # -----------------------------------------------------
  score_translation = function(aa){
    if (nchar(aa) == 0) return(list(stops = Inf, len = 0))
    
    # count '*' stops
    stops = stringr::str_count(aa, fixed("*"))
    # ignore a terminal stop (common in partial CDS)
    if (endsWith(aa, "*")) stops = stops - 1
    
    list(stops = stops, len = nchar(aa))
  }
  
  # -----------------------------------------------------
  # Helper: choose best frame (1/2/3) for one sequence
  #   Criteria:
  #     1) minimum internal stops
  #     2) maximum translated length
  # -----------------------------------------------------
  best_frame = function(dna_str){
    
    # translate all three frames
    a1 = translate_frame(dna_str, 1)
    a2 = translate_frame(dna_str, 2)
    a3 = translate_frame(dna_str, 3)
    
    # score each translation
    s1 = score_translation(a1)
    s2 = score_translation(a2)
    s3 = score_translation(a3)
    
    # pick the best according to (stops asc, length desc)
    df = data.frame(
      frame = c(1,2,3),
      stops = c(s1$stops, s2$stops, s3$stops),
      len   = c(s1$len,   s2$len,   s3$len)
    )
    df = df[order(df$stops, -df$len), ]
    best = df$frame[1]
    
    # build the nucleotide sequence in the best frame (trimmed to codons)
    dna_in_frame = substr(toupper(dna_str), best, nchar(dna_str))
    dna_in_frame = trim_to_full_codons(dna_in_frame)
    
    # AA translation for the chosen frame
    aa_best = list(a1,a2,a3)[[best]]
    
    list(aa = aa_best, dna = dna_in_frame)
  }
  
  # -----------------------------------------------------
  # Helper: back-translate AA alignment to codon alignment
  #   - gaps ('-') in AA alignment become '---' codon gaps
  #   - each AA consumes one codon from dna_in_frame
  # -----------------------------------------------------
  backtranslate = function(aa_aln, dna){
    
    if (nchar(dna) < 3) return("")
    
    # split dna into codons
    starts = seq(1, nchar(dna), by = 3)
    codons = stringr::str_sub(dna, starts, starts+2)
    
    out = character(0)
    k = 1  # codon index
    
    for (ch in strsplit(aa_aln, "")[[1]]){
      
      if (ch == "-"){
        # AA gap -> codon gap
        out = c(out, "---")
      } else {
        # AA residue -> take the next codon
        if (k <= length(codons)){
          out = c(out, codons[k])
          k = k + 1
        } else {
          # alignment longer than sequence (rare edge case): pad with codon gaps
          out = c(out, "---")
        }
      }
    }
    
    paste0(out, collapse="")
  }
  
  # =====================================================
  # 1) Read nucleotide sequences
  # =====================================================
  dna = readDNAStringSet(in_fa)
  dna_strs = as.character(dna)
  names(dna_strs) = names(dna)
  
  # =====================================================
  # 2) Choose best frame per sequence + translate
  # =====================================================
  bf = lapply(dna_strs, best_frame)
  
  # AA strings (for MAFFT) + trimmed in-frame DNA strings (for back-translation)
  aa_strs  = vapply(bf, `[[`, character(1), "aa")
  dna_trim = vapply(bf, `[[`, character(1), "dna")
  
  # Write AA sequences to FASTA
  aa_set = AAStringSet(aa_strs)
  names(aa_set) = names(aa_strs)
  writeXStringSet(aa_set, aa_fa)
  
  # =====================================================
  # 3) Align amino acids with MAFFT
  # =====================================================
  if (!nzchar(Sys.which("mafft")))
    stop("MAFFT not found in PATH.")
  
  message("Running MAFFT...")
  
  # --auto lets MAFFT choose an appropriate strategy; --thread -1 uses all cores
  aln_txt = system2(
    "mafft",
    c("--auto", "--thread", threads, shQuote(aa_fa)),
    stdout = TRUE
  )
  
  # Save MAFFT output as AA alignment FASTA
  writeLines(aln_txt, aa_aln_fa)
  
  # Read AA alignment back in
  aa_aln = readAAStringSet(aa_aln_fa)
  aa_aln_str = as.character(aa_aln)
  
  # =====================================================
  # 4) Back-translate AA alignment to codon alignment
  # =====================================================
  # Keep only names that exist in both alignment and original sequences
  common = intersect(names(aa_aln_str), names(dna_trim))
  
  codon_str = mapply(
    backtranslate,
    aa_aln_str[common],
    dna_trim[common],
    USE.NAMES = TRUE
  )
  
  # Write codon-aligned DNA FASTA
  codon_set = DNAStringSet(codon_str)
  names(codon_set) = names(codon_str)
  writeXStringSet(codon_set, codon_aln)
  
  message("✓ Finished ", locus)
  message("AA alignment:    ", aa_aln_fa)
  message("Codon alignment: ", codon_aln)
  
  # Return file paths invisibly (nice for scripting)
  invisible(list(aa_aln = aa_aln_fa, codon_aln = codon_aln))
}