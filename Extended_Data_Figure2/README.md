# Extended Data Figure 2 — Phylogenetic signal in host use (Fritz & Purvis D)

Code and data for Extended Data Figure 2: phylogenetic signal in host use across the
Danaini butterfly phylogeny, quantified with the Fritz & Purvis (2010) **D-statistic**
(`caper::phylo.d`), the appropriate metric for a binary trait.

A single notebook, built on the same prepared inputs (tree + host associations) as
the Figure 1 ancestral state reconstruction, so the signal test and the
reconstruction rest on the same data:

- **`code/1_hostplant_phylogenetic_signal_Dstatistic.Rmd`** — the D-statistic for
  milkweed-clade host use (the headline test) and for five further host groups (the
  non-milkweed Apocynaceae subfamilies — combined and each alone — and Moraceae),
  with a null-distribution figure and a trait-on-tree figure for each group.

The notebook runs with the working directory at `code/`; all paths are relative to
that folder.

## Directory layout

```
Extended_Data_Figure2/
├── README.md
├── code/
│   └── 1_hostplant_phylogenetic_signal_Dstatistic.Rmd
├── input_data/
│   ├── Danaini_in_HOSTS.newick          # Danaini species tree (from Extended Data Fig. 1)
│   ├── danaini_reference_order.newick   # curated tip order shared with Ext. Fig. 1 / Figure 1
│   └── Ext_fig1_host_associations.csv   # butterfly × plant-genus associations w/ subfamily & family (Ext. Fig. 1)
├── output_data/
│   ├── fritz_purvis_D_milkweed_phylo_signal.txt     # D-stat summary, milkweed clade
│   ├── fritz_purvis_D_<group>_phylo_signal.txt      #   one per further host group
│   └── fritz_purvis_D_host_group_summary.csv        #   D, P(D=1), P(D=0) across all groups
└── visuals/
    ├── hostplant_phylo_signal_D_nulldist[_<group>].png     # D vs. random/Brownian null distributions
    └── hostplant_phylo_signal_trait_tree[_<group>].png     # host-group use mapped on the phylogeny
```

## Inputs (prepared in Extended Data Figure 1)

- **`Danaini_in_HOSTS.newick`** — the filtered Danaini ML tree (56 tips, including
  the outgroups *Mechanitis polymnia* and *Heliconius erato*).
- **`Ext_fig1_host_associations.csv`** — the butterfly-species × plant-genus
  associations shown in Extended Data Fig. 1, each carrying its plant subfamily and
  family; every row is a documented (present) association. Provides the per-group
  presence data for each binary trait.
- **`danaini_reference_order.newick`** — the curated tip order (shared with
  Extended Data Fig. 1 and Figure 1); the trait-on-tree figures are rotated to this
  order so they match the orientation used elsewhere.

A tree tip with no association in a given host group (including the outgroups, and
*Anetia thirza*, which is absent from the association table) is coded `0` for that
group. The D-statistic itself is order-independent; the reference order only sets the
figure orientation.

## Running

```r
setwd("<path>/Extended_Data_Figure2/code")
rmarkdown::render("1_hostplant_phylogenetic_signal_Dstatistic.Rmd")
```

R package dependencies: `tidyverse`, `ape`, `caper`, `ggtree`.
