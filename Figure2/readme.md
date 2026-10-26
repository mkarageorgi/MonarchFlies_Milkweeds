# Figure 2 — Host use, ATPα ancestral states, and ATPase resistance

Code and data for Figure 2. Three notebooks in `code/`, run with the working
directory at `code/` (all paths are relative to that folder):

- **`1_host_use_focused_ASR.Rmd`** — the three-state host-use reconstruction
  (none / milkweed clade / native-range *Asclepias*) restricted to the Danaini
  with a known ATPα genotype, with an inner ATPα strip. The reconstruction is
  **not** refit here: it is read from Figure 1's output, and each focused internal
  node reuses the whole-tree probability of the matching full-tree node (the MRCA
  of the same tips). **Figure 1 must be rendered first** (see below).

- **`2_atpalpha_ASR.Rmd`** — maximum-likelihood ancestral state reconstruction of
  the ATPα functional states across the Danaini, plus the sensitivity analyses.
  Each taxon's state is the three-residue code at ATPα positions 1 / 9 / 12 (VSH,
  VSN, LSN, QAN, LAN). Sections:
  1. **Reconstruction** — IQ-TREE, LG+Γ, marginal (empirical-Bayes), rooted on the
     outgroups *Heliconius erato* and *Mechanitis polymnia*; cross-checked against
     a `phangorn`-LG reconstruction.
  2. **Bootstrap** — the reconstruction repeated across the 1000 ultrafast-boot
     trees of the Danaini backbone (robustness to phylogenetic uncertainty).
  3. **Jackknife** — each sequenced taxon dropped in turn (robustness to taxon
     sampling).
  4. **Coverage** — the reconstruction placed on the full Danaini tree with
     unsampled taxa entered as missing data (where states are supported vs
     interpolated).
  5. **AltAll note** — why the residue-level robustness test is uninformative for
     this dataset.

- **`3_butterfly_vs_fly_IC50s.Rmd`** — ouabain IC₅₀ dose–response analysis
  comparing Na⁺/K⁺-ATPase resistance of milkweed butterflies (Petschenka et al.
  2013) with engineered *Drosophila melanogaster* ATPα knock-in lines (Karageorgi
  et al. 2019). Fits one four-parameter log-logistic curve per biological
  preparation, then compares IC₅₀ values by ATPα state. Produces the supplementary
  IC₅₀ table and the IC₅₀-by-state figure.

## Directory layout

```
Figure2/
├── README.md
├── code/
│   ├── 1_host_use_focused_ASR.Rmd
│   ├── 2_atpalpha_ASR.Rmd
│   └── 3_butterfly_vs_fly_IC50s.Rmd
├── input_data/
│   ├── Danaini_in_HOSTS.newick               # Danaini species tree (from Ext. Data Fig. 1)
│   ├── danaini_reference_order.newick        # curated tip order shared with Ext. Fig. 1 / Figure 1
│   ├── danaini_host_use_states.csv           # per-species host-use states + ATPα (Ext. Fig. 1 notebook 5)
│   ├── alignment.fasta                        # ATPα functional-motif alignment (incl. 4 congener/subspecies-sourced tips)
│   ├── Petschenka2013_NaKATPase_raw_absorbance.csv   # raw colorimetric assay data (Dryad)
│   └── Karageorgi2019_Fig3a_NaKATPase_activity.csv   # Fig. 3a source data (engineered Drosophila)
├── output_data/
│   ├── host_use_focused_ASR_node_probabilities.csv   # (notebook 1) focused-subset node probabilities
│   ├── atpalpha_ASR_node_posteriors.csv              # (notebook 2) IQ-TREE per-node posteriors
│   ├── atpalpha_ASR_bootstrap_robustness.csv         # (notebook 2) per-node state consistency across boot trees
│   ├── atpalpha_ASR_jackknife_stability.csv          # (notebook 2) per-node leave-one-out stability
│   ├── atpalpha_ASR_coverage.csv                     # (notebook 2) full-tree per-node coverage + support
│   ├── atpalpha_ASR_full_tree_posteriors.csv         # (notebook 2) full-tree per-node posteriors
│   ├── iqtree_asr/                                   # (notebook 2) IQ-TREE run files (.state, .treefile, …)
│   └── Supplementary_Table_ATPase_IC50_replicates.xlsx  # (notebook 3) per-preparation IC50 estimates
└── visuals/
    ├── host_use_focused_ASR.png              # (notebook 1) ATPα-known subset reconstruction
    ├── atpalpha_ASR_tree.png                 # (notebook 2) IQ-TREE ancestral states
    ├── atpalpha_ASR_bootstrap.png            # (notebook 2) bootstrap state consistency
    ├── atpalpha_ASR_jackknife.png            # (notebook 2) jackknife state stability
    ├── atpalpha_ASR_coverage.png             # (notebook 2) full-tree coverage diagnostic
    └── IC50_by_ATPase_state_means.png        # (notebook 3) IC50 by ATPα state
```

## External dependencies

**Notebook 1 depends on Figure 1.** It reads the whole-tree reconstruction from
`../../Figure1/output_data/host_use_milkweed_asclepias_ASR_node_probabilities.csv`.
Those probabilities are keyed by the full tree's internal node numbers, so the
tree is rebuilt here with the **same** steps and the **same** `SUPPORT_CUTOFF`
(70) as Figure 1. If the guard in its Section 3 fails, re-render Figure 1, then
this notebook, from the same inputs. The three host-use input files are copied
from `../Extended_Data_Figure1/` so the directory is self-contained apart from the
reconstruction it reads from Figure 1.

**Notebook 2 requires IQ-TREE and the Danaini bootstrap trees.** IQ-TREE is
called via `system2()`; set the `IQTREE_BIN` environment variable to your binary,
otherwise it defaults to `/opt/anaconda3/envs/danaini-phylo/bin/iqtree` and falls
back to `iqtree`/`iqtree2`/`iqtree3` on `PATH`. The bootstrap section reads
`../../Extended_Data_Figure1/output_data/danaini_tree/supermatrix_ge1/supermatrix_ge1.ufboot`.

## Running

```r
# Notebook 1 — render Figure 1 first (produces the reconstruction it reads):
setwd("<path>/Figure1/code")
rmarkdown::render("1_host_use_ancestral_state_reconstruction.Rmd")
setwd("<path>/Figure2/code")
rmarkdown::render("1_host_use_focused_ASR.Rmd")

# Notebook 2 — ATPα ancestral states + sensitivity (needs IQ-TREE on PATH/IQTREE_BIN):
rmarkdown::render("2_atpalpha_ASR.Rmd")

# Notebook 3 — IC50 comparison:
rmarkdown::render("3_butterfly_vs_fly_IC50s.Rmd")
```

## R package dependencies

- **Notebook 1:** `tidyverse`, `ape`, `ggtree`, `ggnewscale`
- **Notebook 2:** `ape`, `phangorn`, `seqinr`, `tidyverse` (plus the IQ-TREE binary)
- **Notebook 3:** `tidyverse`, `drc`, `openxlsx`, `purrr`, `tidyr`
