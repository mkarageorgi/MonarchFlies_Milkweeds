## Figure 1

This code reproduces Figure 1.

## 4. Phylogenetic signal in milkweed-clade host use

4_hostplant_phylogenetic_signal.Rmd estimates phylogenetic signal in milkweed-clade host use across Danaini butterflies using Pagel’s λ. Milkweed-clade use was coded as a binary trait indicating whether each butterfly species has at least one recorded larval host plant in Asclepiadoideae or Secamonoideae.

## 5. Host-use ancestral state reconstruction and Figure 1

5_host_use_ancestral_state_reconstruction.Rmd performs ancestral state reconstruction of Danaini host-plant use and generates Figure 1. The analysis uses the curated host-plant association dataset and the Danaini phylogeny.

For the broad Danaini analysis, host use was coded as a three-state trait: `none`, `milkweed`, or `asclepias`. Species were coded as using the milkweed clade if they had at least one recorded larval host in Asclepiadoideae or Secamonoideae, and as using *Asclepias* if they had at least one recorded larval host in the genus *Asclepias*. Ancestral states were reconstructed using a maximum-likelihood Mk model with equal transition rates.

The R Markdown file also generates a focused subset of species with both ATPα state information and host-use data. For this subset, *Asclepias* use was reconstructed as a binary trait and visualized alongside observed ATPα resistance states.



Code and analysis by Marianna Karageorgi & Naama Weksler.
