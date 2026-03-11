# Taxonomy_clean

# R approach

GBIF_append.R
Flexible R script for harmonising taxonomic data within an arbitrary table of data and appending GBIF results

Usage follows the format:
GBIF_append.R -f input_file -c column_name(s) -out output

taxon_SQL_kingdom.R

Specific script for harmonising the GloBI taxonomy to the GBIF backbone (rather than multiple sources). Uses a mixture of local SQL database and API calls to balance speed and flexibility.

# Python approach

R script seemed fairly slow at the large scales required so there is a faster Python approach.

You can run with GNIF taxon file, input file and the name (header) of the column you want the GBIF matching to run on. Like so:

python GBIF_flexiMatch_prefixer.py ../External_Data/GBIF/Taxon.tsv ./AVONET1_BirdLife.csv.gz Species1 > GBIF_avonet.csv
