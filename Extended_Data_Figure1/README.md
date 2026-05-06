## Extended Data Figure 1

This code reproduces Extended Data Figure 1, which shows host-plant associations for 55 Danaini butterflies and two Danainae outgroups. Larval host records were compiled mainly from Ackery (1988), with additional HOSTS database records. The figure combines the Danaini phylogeny generated here with published Apocynaceae and Moraceae genus-level phylogenies, cardenolide occurrence, and ATPα resistance states.

### 1. Host-plant data curation
Host-plant data used in 1_prepare_hostplant_data.Rmd were manually curated from the sources listed in the Rmd file, primarily Ackery (1988), with additional records from the HOSTS database: https://data.nhm.ac.uk/dataset/hosts. The HOSTS database file, resource.csv, is not included in this repository but can be requested from the Natural History Museum Data Portal.

The 1_prepare_hostplant_data.Rmd generates:

plant_butterfly_long_all.csv: Curate dataset containing butterfly taxonomy, host-plant taxonomy, host-association records, cardenolide occurrence, and ATPα state information.

Code and analysis by Marianna Karageorgi.

