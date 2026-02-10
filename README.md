# Taxonomy_clean

GBIF_append.R
Flexible R script for harmonising taxonomic data within an arbitrary table of data and appending GBIF results

Usage follows the format:
GBIF_append.R -f input_file -c column_name(s) -out output

taxon_SQL_kingdom.R

Specific script for harmonising the GloBI taxonomy to the GBIF backbone (rather than multiple sources). Uses a mixture of local SQL database and API calls to balance speed and flexibility.

