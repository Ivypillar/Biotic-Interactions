library(GIFT) 
library(data.table) 
library(ggplot2) 
library(rgbif) 
#library(traitdataform) 
library(purrr) 
library(optparse)
library(dplyr)
library(taxize)
library(taxizedb)
library(pbapply)
library(retry)
library(sf) 
library(dplyr)
#library(fuzzyjoin)

# setting up taxonomy resources locally
#interactions <- fread("~/Desktop/PostDoc_Ecology/Interactions.csv.gz")
#interactions <- fread("../External_Data/globi_IntDataProd_0.8/interactions.csv.gz")

db_download_gbif(verbose = TRUE, overwrite = FALSE)
src_gbif <- src_gbif()
#gbif <- src_gbif %>% tbl("gbif")
#options(dplyr.width = Inf)
#sql_collect(src_gbif, "select * from gbif limit 5")

GBIF_SQL_lookup <- function(input_temp){
  input<-data.frame(t(input_temp))
  SQL_source_string <- paste0("select * from gbif where canonicalName = '",gsub('[[:punct:] ]+',' ',input$sourceTaxonName),"'")
  SQL_target_string <- paste0("select * from gbif where canonicalName = '",gsub('[[:punct:] ]+',' ',input$targetTaxonName),"'")
  #print(paste(input$sourceTaxonName,"<->",input$targetTaxonName))
  ## Try local SQL first for speed.
  source_collect <- subset(sql_collect(src_gbif, SQL_source_string),taxonomicStatus=="accepted")
  #print(source_collect)
  source_collect <-  tibble(GBIF_acceptedID=ifelse("taxonID" %in% colnames(source_collect),source_collect[,"taxonID"],NA),
                            GBIF_taxonRank=ifelse("taxonRank" %in% colnames(source_collect),source_collect[,"taxonRank"],NA),
                            GBIF_kingdom=ifelse("kingdom" %in% colnames(source_collect),source_collect[,"kingdom"],NA),
                            GBIF_phylum=ifelse("phylum" %in% colnames(source_collect),source_collect[,"phylum"],NA),
                            GBIF_class=ifelse("class" %in% colnames(source_collect),source_collect[,"class"],NA),
                            GBIF_order=ifelse("order" %in% colnames(source_collect),source_collect[,"order"],NA),
                            GBIF_family=ifelse("family" %in% colnames(source_collect),source_collect[,"family"],NA),
                            GBIF_genus=ifelse("genus" %in% colnames(source_collect),source_collect[,"genus"],NA),
                            GBIF_species=ifelse("specificEpithet" %in% colnames(source_collect),source_collect[,"specificEpithet"],NA)
  )
  
  target_collect <- subset(sql_collect(src_gbif, SQL_target_string),taxonomicStatus=="accepted")
  target_collect <-  tibble(GBIF_acceptedID=ifelse("taxonID" %in% colnames(target_collect),target_collect[,"taxonID"],NA),
                            GBIF_taxonRank=ifelse("taxonRank" %in% colnames(target_collect),target_collect[,"taxonRank"],NA),
                            GBIF_kingdom=ifelse("kingdom" %in% colnames(target_collect),target_collect[,"kingdom"],NA),
                            GBIF_phylum=ifelse("phylum" %in% colnames(target_collect),target_collect[,"phylum"],NA),
                            GBIF_class=ifelse("class" %in% colnames(target_collect),target_collect[,"class"],NA),
                            GBIF_order=ifelse("order" %in% colnames(target_collect),target_collect[,"order"],NA),
                            GBIF_family=ifelse("family" %in% colnames(target_collect),target_collect[,"family"],NA),
                            GBIF_genus=ifelse("genus" %in% colnames(target_collect),target_collect[,"genus"],NA),
                            GBIF_species=ifelse("specificEpithet" %in% colnames(target_collect),target_collect[,"specificEpithet"],NA)
  )
  #print(unlist(source_collect))
  
  # If SQL doesn't give a clean match, try the API fuzzy match and constrain by confidence and Order
  if(length(as.integer(unlist(source_collect["GBIF_acceptedID"])))!=1){
    #print("Source: No SQL entry, trying API")
    #print(input)
    
    source_api_call_init <- get_gbifid_(input$sourceTaxonName,messages=F)[[1]]
    
    
    
    #print(paste(input$sourceTaxonKingdomName,input$sourceTaxonOrderName,input$sourceTaxonName))
    #print(source_api_call_init)
    if(length(source_api_call_init)!=0 & "kingdom" %in% colnames(source_api_call_init)){
      
      if(input$sourceTaxonKingdomName=="" | is.na(input$sourceTaxonKingdomName) ){
        #print("if")
        source_api_call <- source_api_call_init[1,]
        
      }
      else if(length(subset(source_api_call_init,kingdom==input$sourceTaxonKingdomName)[,1])==0  & input$sourceTaxonRank=="species" & input$sourceTaxonGenusName!=""){
        #print("if")
        source_api_call_init <- get_gbifid_(input$sourceTaxonGenusName,messages=F)[[1]]
        if((length(subset(source_api_call_init,kingdom==input$sourceTaxonKingdomName)[,1])==0 )){
          return(NULL)
        }
        else{
          source_api_call <- subset(source_api_call_init,(kingdom==input$sourceTaxonKingdomName ) )[1,]
        }
      }
      else if((length(subset(source_api_call_init,kingdom==input$sourceTaxonKingdomName)[,1])==0 ) & input$sourceTaxonRank!="species"){
        return(NULL)
      }
      else if(source_api_call_init$kingdom[1]!=input$sourceTaxonKingdomName ){
        return(NULL)
      } 
      else{
        #print("else")
        
        source_api_call <- subset(source_api_call_init,(kingdom==input$sourceTaxonKingdomName ) )[1,]
        #print(source_api_call)
      }
      
      #print(source_api_call)
      if(source_api_call$status=="ACCEPTED"){
        source_collect <-  tibble(GBIF_acceptedID=source_api_call$usagekey,
                                  GBIF_taxonRank=source_api_call$rank,
                                  GBIF_kingdom=ifelse("kingdom" %in% colnames(source_api_call),source_api_call$kingdom,NA),
                                  GBIF_phylum=ifelse("phylum" %in% colnames(source_api_call),source_api_call$phylum,NA),
                                  GBIF_class=ifelse("class" %in% colnames(source_api_call),source_api_call$class,NA),
                                  GBIF_order=ifelse("order" %in% colnames(source_api_call),source_api_call$order,NA),
                                  GBIF_family=ifelse("family" %in% colnames(source_api_call),source_api_call$family,NA),
                                  GBIF_genus=ifelse("genus" %in% colnames(source_api_call),source_api_call$genus,NA),
                                  GBIF_species=ifelse("species" %in% colnames(source_api_call),source_api_call$species,NA)
        )
      }
      #source_collect <- tibble(taxonID=NA,taxonRank=NA)
      if(source_api_call$status=="SYNONYM"){
        source_collect <-  tibble(GBIF_acceptedID=source_api_call$acceptedusagekey,
                                  GBIF_taxonRank=source_api_call$rank,
                                  GBIF_kingdom=ifelse("kingdom" %in% colnames(source_api_call),source_api_call$kingdom,NA),
                                  GBIF_phylum=ifelse("phylum" %in% colnames(source_api_call),source_api_call$phylum,NA),
                                  GBIF_class=ifelse("class" %in% colnames(source_api_call),source_api_call$class,NA),
                                  GBIF_order=ifelse("order" %in% colnames(source_api_call),source_api_call$order,NA),
                                  GBIF_family=ifelse("family" %in% colnames(source_api_call),source_api_call$family,NA),
                                  GBIF_genus=ifelse("genus" %in% colnames(source_api_call),source_api_call$genus,NA),
                                  GBIF_species=ifelse("species" %in% colnames(source_api_call),source_api_call$species,NA))
      }
    }
  }
  #print(length(as.integer(unlist(target_collect["GBIF_acceptedID"]))))
  if(length(as.integer(unlist(target_collect["GBIF_acceptedID"])))!=1){
    #print("Target: No SQL entry, trying API")
    #print(input)
    
    target_api_call_init <- get_gbifid_(input$targetTaxonName,messages=F)[[1]]
    
    #check if gbif higher taxonomy matches GloBI data, if not, I have found that the species doesn't exist but the correct genus is returned  if that is used solely.
    
    
    #target_api_call_init <- get_gbifid_(input$targetTaxonName,messages=F)[[1]]
    #print(target_api_call_init)
    #print(paste(input$targetTaxonKingdomName[1],input$targetTaxonPhylumName,input$targetTaxonName[1]))
    #print(paste("kingdom=",target_api_call_init$kingdom[1],"targetTaxonKingdomName=",input$targetTaxonKingdomName))
    
    #print(subset(target_api_call_init,kingdom==input$targetTaxonKingdomName))
    #check that the api call has returned a result.
    if(length(target_api_call_init)!=0 & "kingdom" %in% colnames(target_api_call_init)){
      
      #print(input$targetTaxonKingdomName)
      #if there's no higher GloBI taxon info then we do a general search
      if(input$targetTaxonKingdomName=="" | is.na(input$targetTaxonKingdomName) ){
        #print("if")
        target_api_call <- target_api_call_init[1,]
        
      }
      else if((length(subset(target_api_call_init,kingdom==input$targetTaxonKingdomName)[,1])==0 ) & input$targetTaxonRank=="species" & input$targetTaxonGenusName!=""){
        #print("elseif1")
        #print(input$targetTaxonGenusName)
        target_api_call_init <- get_gbifid_(input$targetTaxonGenusName,messages=F)[[1]]
        #print(target_api_call_init)
        if((length(subset(target_api_call_init,kingdom==input$targetTaxonKingdomName)[,1])==0 )){
          return(NULL)
        }
        else{
          target_api_call <- subset(target_api_call_init,(kingdom==input$targetTaxonKingdomName ) )[1,]
        }
      }
      else if((length(subset(target_api_call_init,kingdom==input$targetTaxonKingdomName)[,1])==0 ) & input$targetTaxonRank!="species"){
        #print("elseif2")  
        return(NULL)
      }
      
      else if(target_api_call_init$kingdom[1]!=input$targetTaxonKingdomName){
        #print("elseif3")  
        return(NULL)
        
      } 
      
      #else we make sure that the subset includes the higher taxonomy to avoid misclassification of species with the same/similar names from diferent phyla
      else{
        #print("else")
        target_api_call <- subset(target_api_call_init,(kingdom==input$targetTaxonKingdomName  ) )[1,]
      }
      #if(target_api_call$status=="ACCEPTED"){print("test")}
      #print(target_api_call)
      #print(paste("target",input$targetTaxonKingdomName,input$targetTaxonName,input$sourceTaxonName,input$interactionTypeName))
      if(target_api_call$status=="ACCEPTED"){
        target_collect <-  tibble(GBIF_acceptedID=target_api_call$usagekey,
                                  GBIF_taxonRank=target_api_call$rank,
                                  GBIF_kingdom=ifelse("kingdom" %in% colnames(target_api_call),target_api_call$kingdom,NA),
                                  GBIF_phylum=ifelse("phylum" %in% colnames(target_api_call),target_api_call$phylum,NA),
                                  GBIF_class=ifelse("class" %in% colnames(target_api_call),target_api_call$class,NA),
                                  GBIF_order=ifelse("order" %in% colnames(target_api_call),target_api_call$order,NA),
                                  GBIF_family=ifelse("family" %in% colnames(target_api_call),target_api_call$family,NA),
                                  GBIF_genus=ifelse("genus" %in% colnames(target_api_call),target_api_call$genus,NA),
                                  GBIF_species=ifelse("species" %in% colnames(target_api_call),target_api_call$species,NA))
      }
      #source_collect <- tibble(taxonID=NA,taxonRank=NA)
      if(target_api_call$status=="SYNONYM"){
        target_collect <-  tibble(GBIF_acceptedID=target_api_call$acceptedusagekey,
                                  GBIF_taxonRank=target_api_call$rank,
                                  GBIF_kingdom=ifelse("kingdom" %in% colnames(target_api_call),target_api_call$kingdom,NA),
                                  GBIF_phylum=ifelse("phylum" %in% colnames(target_api_call),target_api_call$phylum,NA),
                                  GBIF_class=ifelse("class" %in% colnames(target_api_call),target_api_call$class,NA),
                                  GBIF_order=ifelse("order" %in% colnames(target_api_call),target_api_call$order,NA),
                                  GBIF_family=ifelse("family" %in% colnames(target_api_call),target_api_call$family,NA),
                                  GBIF_genus=ifelse("genus" %in% colnames(target_api_call),target_api_call$genus,NA),
                                  GBIF_species=ifelse("species" %in% colnames(target_api_call),target_api_call$species,NA))
      }
      #target_collect <- tibble(taxonID=NA,taxonRank=NA)
    }
  }
  #print(paste(source_collect))
  colnames(source_collect) <- paste("source", colnames(source_collect), sep = "_")
  colnames(target_collect) <- paste("target", colnames(target_collect), sep = "_")
  out <- cbind(input,source_collect,target_collect)
  #print(out)
  if(is.na(out$source_GBIF_phylum)==F & is.na(out$target_GBIF_phylum)==F){
  if(out$source_GBIF_phylum=="Tracheophyta" | out$target_GBIF_phylum=="Tracheophyta"){
    return(out)
  }
}
}


#GBIF_harm <-rbindlist(pbapply(problem_record,1, function(x) GBIF_SQL_lookup(x) ))

GBIF_harm <-rbindlist(pbapply(interactions,1, function(x) GBIF_SQL_lookup(x) ))


#
#print(tail(GBIF_harm))


#subset(target_api_call_init, phylum == input$targetTaxonPhylumName)


