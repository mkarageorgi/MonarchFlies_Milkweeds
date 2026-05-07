## Extended Data Figure 1

This code reproduces Extended Data Figure 1, which shows host-plant associations for 55 Danaini butterflies and two Danainae outgroups. Larval host records were compiled mainly from Ackery (1988), with additional HOSTS database records. The figure combines the Danaini phylogeny generated here with published Apocynaceae and Moraceae genus-level phylogenies, cardenolide occurrence, and ATPα resistance states.

### 1. Host-plant data curation
1_prepare_hostplant_data.Rmd prepares the host-plant association dataset. Host-plant records were manually curated from the sources listed in the Rmd file, primarily Ackery (1988), with additional records from the HOSTS database: https://data.nhm.ac.uk/dataset/hosts. The HOSTS database file, resource.csv, is not included in this repository but can be requested from the Natural History Museum Data Portal.

The script generates:

plant_butterfly_long_all.csv, a curated dataset containing butterfly taxonomy, host-plant taxonomy, host-association records, cardenolide occurrence, and ATPα state information.

### 2. Danaini phylogeny reconstruction
2_prepare_danaini_phylogeny.Rmd generates the Danaini phylogeny used in Extended Data Figure 1 & Figure 1. Related scripts are within the scripts folder. The analysis follows the approach of Chazot et al. (2021) and uses available mitochondrial and nuclear gene sequences for Danaini butterflies and two Danainae outgroups.

Sequences were aligned, concatenated, and used to infer a maximum-likelihood phylogeny. The resulting tree was pruned to the 55 Danaini species and two outgroups with host-plant data for subsequent analyses.

### 3. Host-plant association matrix (Extended Data Figure 1)

3_reconstruct_hostplant_associations.Rmd generates Extended Data Figure 1 by combining the Danaini butterfly phylogeny reconstructed in this repository, genus-level Apocynaceae and Moraceae phylogenies, and the curated butterfly–host plant association matrix.

The script first prunes and roots the Danaini phylogeny to include the 55 Danaini species and two Danainae outgroups with host-plant data. For the plant phylogenies, the script uses the time-calibrated Apocynaceae phylogeny from Fishbein et al. (2018), which was pruned to one representative species per genus and then renamed to generate a genus-level Apocynaceae tree. Moraceae relationships were based on a genus-level tree obtained from TimeTree, with genus names changed to current taxonomy where needed.

The Apocynaceae and Moraceae trees were then pruned to the genera represented in the host-plant matrix and grafted into a combined plant phylogeny. The script also generates a butterfly species × plant genus presence/absence matrix and visualizes host use together with plant cardenolide occurrence and butterfly ATPα resistance states.

Code and analysis by Marianna Karageorgi.

