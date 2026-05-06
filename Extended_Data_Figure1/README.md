## Extended Data Figure 1

This code reproduces Extended Data Figure 1, which shows host-plant associations for 55 Danaini butterflies and two Danainae outgroups. Larval host records were compiled mainly from Ackery (1988), with additional HOSTS database records. The figure combines the Danaini phylogeny generated here with published Apocynaceae and Moraceae genus-level phylogenies, cardenolide occurrence, and ATPα resistance states.

### 1. Host-plant data curation
1_prepare_hostplant_data.Rmd prepares the host-plant association dataset. Host-plant records were manually curated from the sources listed in the Rmd file, primarily Ackery (1988), with additional records from the HOSTS database: https://data.nhm.ac.uk/dataset/hosts. The HOSTS database file, resource.csv, is not included in this repository but can be requested from the Natural History Museum Data Portal.

The script generates:

plant_butterfly_long_all.csv, a curated dataset containing butterfly taxonomy, host-plant taxonomy, host-association records, cardenolide occurrence, and ATPα state information.

### 2. Danaini phylogeny reconstruction
2_prepare_danaini_phylogeny.Rmd generates the Danaini phylogeny used in Extended Data Figure 1 & Figure 1. Related scripts are within the scripts folder. The analysis follows the approach of Chazot et al. (2021) and uses available mitochondrial and nuclear gene sequences for Danaini butterflies and two Danainae outgroups.

Sequences were aligned, concatenated, and used to infer a maximum-likelihood phylogeny. The resulting tree was pruned to the 55 Danaini species and two outgroups with host-plant data for subsequent analyses.

Code and analysis by Marianna Karageorgi.

