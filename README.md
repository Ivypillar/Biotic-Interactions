# Taxonomy_clean

# R approach

GBIF_append.R
Flexible R script for harmonising taxonomic data within an arbitrary table of data and appending GBIF results

Usage follows the format:
GBIF_append.R -f input_file -c column_name(s) -out output

taxon_SQL_kingdom.R

Specific script for harmonising the GloBI taxonomy to the GBIF backbone (rather than multiple sources). Uses a mixture of local SQL database and API calls to balance speed and flexibility.

# Python approach

R script seemed fairly slow at the large scales required so there is a faster Python approach. I've included a yaml file of the conda environment, to pre-empt dependency issues:
conda env create -f GBIF_harmonise.yml

You can run with GBIF taxon file, input file and the name (header) of the column you want the GBIF matching to run on, aswell as a secondary backup name column (if this exists) and also supporting phylogeny information to resolve multiple name matches. Like so:

python GBIF_flexiMatch_prefixer.py [GBIF_taxon_file] [Input_file] [Input_taxon_column] [Backup_input_taxon_column] [comma-separated phylogeny columns]  > [output.csv]

python GBIF_flexiMatch_prefixer.py ../External_Data/GBIF/Taxon.tsv ./AVONET1_BirdLife.csv.gz Species1 Species1 Order1,Family1 > GBIF_avonet.csv

python GBIF_flexiMatch_prefixer_author.py ../External_Data/GBIF/Taxon2.tsv ./Pollination_Dispersal_data.csv.gz work_species work_species work_species,work_species work_author > GBIF_GIFT_Pollination_Dispersal_data.csv


