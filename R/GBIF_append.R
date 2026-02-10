# Appends GBIF taxonomy data to specified columns 

# Data from: GBIF 
#  
#  

# Authors: Ivy Whiteford 
# Version: 12-01-2026 

# Packages 
library(GIFT) 
library(data.table) 
library(ggplot2) 
library(sf) 
library(rgbif) 
library(traitdataform) 
library(purrr) 
library(optparse)
# Parameters in the form "script.R input_file.csv column1 column2 …" 


option_list = list(
  make_option(c("-f", "--file"), type="character", default="~/BIDENS/External_Data/globi_IntDataProd_0.8/interactions_loc3.csv", 
              help="dataset file name", metavar="file_string"),
  make_option(c("-c", "--cols"), type="character", default="sourceTaxonSpeciesName targetTaxonSpeciesName", 
              help="Column names to reassign separated by spaces", metavar="column_string"),
  make_option(c("-g", "--globi"), type="logical", default="", 
              help="Resolve missing taxanomic levels in GloBI", metavar="globi_mode"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name", metavar="out_string")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

columns <- strsplit(opt$cols, " ")[[1]]

# test if there is at least one argument: if not, return an error 

input_file <- fread(opt$file) 

#=============================================================================== 


input_file_na  <- input_file %>% mutate_if(is.character, ~na_if(., '')) 
#going to use a for loop here as there shouldn't be many columns.
out_df <- as.data.frame(input_file_na)
if (opt$globi == T) {
  
    
    input_file_na<- input_file_na %>% mutate(source_phylo = coalesce(sourceTaxonSpeciesName, sourceTaxonGenusName, sourceTaxonFamilyName,sourceTaxonOrderName,sourceTaxonClassName,sourceTaxonPhylumName, sourceTaxonKingdomName))
    input_file_na<- input_file_na %>% mutate(target_phylo = coalesce(targetTaxonSpeciesName, targetTaxonGenusName, targetTaxonFamilyName, targetTaxonOrderName, targetTaxonClassName, targetTaxonPhylumName, targetTaxonKingdomName))
  
    print(paste0("GBIF Assigning Column: source_phylo_GBIF"))
    source_phylo_gbif_result <- get_gbif_taxonomy(input_file_na[,source_phylo],subspecies=F,higherrank=T) 
    print(paste0("GBIF Assigning Column: target_phylo_GBIF"))
    target_phylo_gbif_result <- get_gbif_taxonomy(input_file_na[,target_phylo],subspecies=F,higherrank=T) 
    
    colnames(source_phylo_gbif_result) <- paste("source_phylo_GBIF", colnames(source_phylo_gbif_result), sep = "_")
    colnames(target_phylo_gbif_result) <- paste("target_phylo_GBIF", colnames(target_phylo_gbif_result), sep = "_")
    
    out_df<-cbind(out_df,source_phylo_gbif_result,target_phylo_gbif_result)
    
} else 
  for (x in 1:length(columns)) {
    print(paste0("GBIF Assigning Column: ",columns[x]))
    #mapping to species level or higher
    
    gbif_result <- get_gbif_taxonomy(input_file[,eval(as.name(columns[x]))],subspecies=F,higherrank=T) 
    colnames(gbif_result) <- paste(columns[x],"GBIF", colnames(gbif_result), sep = "_")
    out_df<-cbind(out_df,gbif_result)
  }

write.csv(out_df, opt$out, row.names = FALSE)
