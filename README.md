# Zostera-capensis-genomics
Pool-seq data and code for Zostera capensis_South Africa.

The following files contain measures of genomic diversity for Zostera capensis generated in 
Phair et al. 2019 (Anthropogenic pressures negatively impact genomic diversity of the widespread but vulnerable seagrass Zostera capensis):
AR_env_loci_edit.csv - Allelic richness.
he_env_loci_edit1.csv - Expected heterozygosity.
nucleotide_div_env_edit.csv - Nucleotide diversity.

diversity_indices_status.csv contains genomic diversity averaged for each estuary, anthropogentic pressures and the area of submerged macrophytes (from the NBA 2012 - Van niekerk 2012).

env_status.csv contains anthropogenic pressures and ecological category for each estuary (from the NBA 2012 - Van niekerk 2012).

glm.R contains the script used to examine the relationship between the anthropogenic pressures and the measures of genomic diversity. This was carried out using a GLM and post hoc Tukey HSD test.
