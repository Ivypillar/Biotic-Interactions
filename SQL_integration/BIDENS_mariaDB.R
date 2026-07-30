#source ~/.bashrc
#micromamba activate BIDENS_V1

library(RMariaDB)
library(DBI)
library(data.table)
library(dplyr)
library(dbplyr)
library(sf)
library(ggplot2)


con <- dbConnect(RMariaDB::MariaDB(),
                 host = "127.0.0.1",
                 user = "root",
                 password = "new_password",dbname="BIDENS_V1")

#CREATE DATABASE BIDENS_V1;
#USE BIDENS_V1;

#only need to run once to build the mariaDB database!
glonaf_list <- fread("../External_Data/GloNAF/glonaf_list.csv")
glonaf_region <- fread("../External_Data/GloNAF/glonaf_region.csv")
glonaf_taxon_wcvp <- fread("../External_Data/GloNAF/glonaf_taxon_wcvp.csv")
glonaf_flora2 <- fread("../External_Data/GloNAF/glonaf_flora2.csv")
glonaf_TxR <- fread("../External_Data/GloNAF/glonaf_TxR.csv")
GBIF_glonaf_taxon_wcvp_author <- fread("./GBIF_glonaf_taxon_wcvp_author.csv")

#dbWriteTable(con, "glonaf_region", glonaf_region)
#dbWriteTable(con, "glonaf_list", glonaf_list)
#dbWriteTable(con, "glonaf_taxon_wcvp", glonaf_taxon_wcvp)
#dbWriteTable(con, "glonaf_flora2", glonaf_flora2)
#dbWriteTable(con, "glonaf_TxR", glonaf_TxR)
#dbWriteTable(con, "GBIF_glonaf_taxon_wcvp_author", GBIF_glonaf_taxon_wcvp_author)


# Code for loading data from the mariaDB. Reference the tables via the DBI connection
#glonaf_flora2_tbl <- tbl(con, "glonaf_flora2")
#glonaf_list_tbl <- tbl(con, "glonaf_list")
#glonaf_region_tbl <- tbl(con, "glonaf_region")
#glonaf_TxR_tbl <- tbl(con, "glonaf_TxR")
#GBIF_glonaf_taxon_wcvp_author_tbl <- tbl(con, "GBIF_glonaf_taxon_wcvp_author")

my_sf <- read_sf("../External_Data/GloNAF/glonaf_2024_nol/glonaf_2024_nol.shp")

#ggplot(my_sf) + geom_sf(fill = "#69b3a2", color = "white") + theme_void()


#USE BIDENS_V1;

# Link and filter the tables using standard tidyverse syntax
#linked_tbl <- glonaf_list_tbl %>% left_join(glonaf_region_tbl, by = c("region_id" = "id")) 

linked_tbl <- glonaf_flora2_tbl %>% left_join(GBIF_glonaf_taxon_wcvp_author_tbl, by = c("taxon_wcvp_id"="taxon_orig_id") ) 
linked_tbl2 <- linked_tbl %>% left_join(glonaf_region_tbl, by = c("region_id"="id") ) 
linked_data2_coll <- collect(linked_tbl2)
linked_tbl3 <- linked_data2_coll %>% left_join(my_sf, by ="OBJIDsic") 
linked_tbl3$GBIF_species <- paste0(linked_tbl3$GBIF_genericName," ",linked_tbl3$GBIF_specificEpithet)

test_top100 <- fread("Top100_BioticInteractions_withGloNAF_GIFT_family.csv")
Pinus_radiata <-linked_tbl3[linked_tbl3$GBIF_species=="Pinus radiata",]



# Pull the linked data into R
#linked_data <- collect(linked_tbl)


library(readr)
library(tidyr)
library(dplyr)
library(pROC)
library(future)
library(future.apply)
library(kohonen)
library(missSOM)
library(ggrepel)
library(data.table)
library(dplyr)
library(ggalluvial)
#library(GIFT)
#library(reshape2)
library(ggpubr)
library(sf)
library(ggforce)
library(V.PhyloMaker2)

#harmonised <- fread("GBIF_harmonised_GloBI_fast.tsv")
harmonised <- fread("../harmonised/GBIF_harmonised_GloBI_fast.tsv")
harmonised$sourceName <- ifelse(harmonised$GBIF_source_specificEpithet=="",paste0(harmonised$GBIF_source_genericName," sp"),paste0(harmonised$GBIF_source_genericName," ",harmonised$GBIF_source_specificEpithet))
harmonised$targetName <- ifelse(harmonised$GBIF_target_specificEpithet=="",paste0(harmonised$GBIF_target_genericName," sp"),paste0(harmonised$GBIF_target_genericName," ",harmonised$GBIF_target_specificEpithet))


interaction_simplification_table <- data.frame(
  globi_term=c("eats",
               "eatenBy",
               "preysOn",
               "preyedUponBy",
               "kills",
               "killedBy",
               "parasiteOf",
               "hasParasite",
               "endoparasiteOf",
               "hasEndoparasite",
               "ectoparasiteOf",
               "hasEctoparasite",
               "parasitoidOf",
               "hasParasitoid",
               "hostOf",
               "hasHost",
               "pollinates",
               "pollinatedBy",
               "pathogenOf",
               "allelopathOf",
               "hasPathogen",
               "vectorOf",
               "hasVector",
               "dispersalVectorOf",
               "hasDispersalVector",
               "rootparasiteOf",
               "hemiparasiteOf",
               "hasHabitat",
               "createsHabitatFor",
               "epiphyteOf",
               "hasEpiphyte",
               "providesNutrientsFor",
               "acquiresNutrientsFrom",
               "symbiontOf",
               "mutualistOf",
               "commensalistOf",
               "flowersVisitedBy",
               "visitsFlowersOf",
               "ecologicallyRelatedTo",
               "coOccursWith",
               "coRoostsWith",
               "interactsWith",
               "adjacentTo"),
  simplified_term =
    c("Herbivory",
      "Herbivory",
      "Unclear",
      "Unclear",
      "Unclear",
      "Unclear",
      "Parasitism",
      "Parasitism",
      "Parasitism",
      "Parasitism",
      "Parasitism",
      "Parasitism",
      "Parasitism",
      "Parasitism",
      "Host",
      "Host",
      "Pollination",
      "Pollination",
      "Pathogen",
      "Allelopath",
      "Pathogen",
      "Pathogen",
      "Pathogen",
      "Dispersal",
      "Dispersal",
      "Parasitism",
      "Parasitism",
      "Niche_Creation",
      "Niche_Creation",
      "Commensalism",
      "Commensalism",
      "Symbiosis",
      "Symbiosis",
      "Symbiosis",
      "Mutualism",
      "Commensalism",
      "Unclear",
      "Unclear",
      "Unclear",
      "Co-location",
      "Unclear",
      "Unclear",
      "Co-location"))

#df1 %>% 


harmonised_int <- harmonised %>% 
  inner_join(interaction_simplification_table, by = c("interactionTypeName" = "globi_term")) %>%
  mutate(interaction_simple = simplified_term)  %>%
  filter(if_any(all_of(c("GBIF_source_phylum", "GBIF_target_phylum")), ~grepl("Tracheophyta", .)))


harmonised_comp <- subset(harmonised_int,GBIF_source_phylum!="" & GBIF_target_phylum!="" & GBIF_source_phylum!="NA" & GBIF_target_phylum!="NA")
#harmonised_comp$order_pair <- paste0(harmonised_comp$GBIF_source_order,"-",harmonised_comp$GBIF_target_order)
#harmonised_comp$phylum_pair <- paste0(harmonised_comp$GBIF_source_phylum,"-",harmonised_comp$GBIF_target_phylum)
#harmonised_comp$class_pair <- paste0(harmonised_comp$GBIF_source_class,"-",harmonised_comp$GBIF_target_class)

harmonised_comp_rearr <- harmonised_comp %>% 
  mutate(new_target_kingdom=ifelse(GBIF_target_phylum=="Tracheophyta",GBIF_source_kingdom,GBIF_target_kingdom),
         new_source_kingdom=ifelse(GBIF_target_phylum=="Tracheophyta",GBIF_target_kingdom,GBIF_source_kingdom)) %>%
  mutate(new_target_phylum=ifelse(GBIF_target_phylum=="Tracheophyta",GBIF_source_phylum,GBIF_target_phylum),
         new_source_phylum=ifelse(GBIF_target_phylum=="Tracheophyta",GBIF_target_phylum,GBIF_source_phylum)) %>% 
  mutate(new_target_class=ifelse(GBIF_target_phylum=="Tracheophyta",GBIF_source_class,GBIF_target_class),
         new_source_class=ifelse(GBIF_target_phylum=="Tracheophyta",GBIF_target_class,GBIF_source_class)) %>% 
  mutate(new_target_order=ifelse(GBIF_target_phylum=="Tracheophyta",GBIF_source_order,GBIF_target_order),
         new_source_order=ifelse(GBIF_target_phylum=="Tracheophyta",GBIF_target_order,GBIF_source_order)) %>% 
  mutate(new_target_family=ifelse(GBIF_target_phylum=="Tracheophyta",GBIF_source_family,GBIF_target_family),
         new_source_family=ifelse(GBIF_target_phylum=="Tracheophyta",GBIF_target_family,GBIF_source_family)) %>% 
  mutate(new_target_genericName=ifelse(GBIF_target_phylum=="Tracheophyta",GBIF_source_genericName,GBIF_target_genericName),
         new_source_genericName=ifelse(GBIF_target_phylum=="Tracheophyta",GBIF_target_genericName,GBIF_source_genericName)) %>% 
  mutate(new_target_specificEpithet=ifelse(GBIF_target_phylum=="Tracheophyta",GBIF_source_specificEpithet,GBIF_target_specificEpithet),
         new_source_specificEpithet=ifelse(GBIF_target_phylum=="Tracheophyta",GBIF_target_specificEpithet,GBIF_source_specificEpithet)) %>% 
  mutate(new_targetName=ifelse(GBIF_target_phylum=="Tracheophyta",sourceName,targetName),
         new_sourceName=ifelse(GBIF_target_phylum=="Tracheophyta",targetName,sourceName)) %>% 
  mutate(new_target_taxonID=ifelse(GBIF_target_phylum=="Tracheophyta",GBIF_source_taxonID,GBIF_target_taxonID),
         new_source_taxonID=ifelse(GBIF_target_phylum=="Tracheophyta",GBIF_target_taxonID,GBIF_source_taxonID))




#head(harmonised_comp_rearr[,c(116,118,120,122,124,126,128,117,119,121,123,125,127,129,113)])

GIDIAS <- fread("../External_Data/GIDIAS/GIDIAS_20250417_machine_read.csv",sep=",")
GIDIAS_mag <- GIDIAS[,c("Verified.Name.GBIF.Taxon","magnitude.Nature","direction.Nature","Country.Location")]
GIDIAS_max <- aggregate(magnitude.Nature ~ Verified.Name.GBIF.Taxon , data = GIDIAS_mag, FUN = function(x) c(mx = max(x), cnt = length(x), sum = sum(x) ))

#avian trait data
avonet <- fread("GBIF_avonet.csv")
avonet$GBIF_taxonID <- as.character(avonet$GBIF_taxonID)

#Gift plant trait data
GBIF_GIFT_Pollination <- fread("GBIF_GIFT_Pollination_Dispersal_data.csv")
GBIF_GIFT_Pollination$GBIF_taxonID <- as.character(GBIF_GIFT_Pollination$GBIF_taxonID)

#GloNAF data
GBIF_GloNAF <- fread("GBIF_glonaf_taxon_wcvp_author.csv")
GBIF_GloNAF_uniq <- unique(GBIF_GloNAF$GBIF_taxonID)


#TRY database
#28	Dispersal syndrome	477217	18582	457440	52887
#231	Dispersal unit type	196259
#207	Flower color	20434	202	20434	10587
#3117	Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or exclu	249181	206176	228304	16599
#677	Leaf emergences (pubescence, pruinescence, hairs, spines, thorns)	33485	28055	33195	3412
#935	Leaf herbivore species	1680		1680	60
#37	Leaf phenology type	240394	128205	207819	29770
#7	Mycorrhiza type	118390	57195	117528	8028
#22	Photosynthesis pathway	146169
#42	Plant growth form	2371038	625152	2333458	233812
#3400	Plant growth form simple consolidated	213372
#77	Plant growth rate relative (plant relative growth rate, RGR)	12682
#3107	Plant height generative	41669	20363	39431	3734
#3106	Plant height vegetative	368479	268607	340236	35466
#343	Plant life form (Raunkiaer life form)	117350	63988	60226	16558
#335	Plant reproductive phenology timing (flowering time)	100818	51074	100578	10856
#38	Plant woodiness	239889	90316	211060	
#29	Pollination syndrome	33745	1457	32625	15879
#3096	Species habitat characterization: vegetation type	85202		85202	12809






#full_df = base::merge(harmonised_comp_rearr, GBIF_GIFT_Pollination, by.x = c('new_source_taxonID'), by.y= c('GBIF_taxonID'))

##Plant_merges

#harmonised_comp_rearr %>% left_join(GBIF_GloNAF, by = join_by(new_source_taxonID == GBIF_taxonID))

#full_df = base::merge(harmonised_comp_rearr, GBIF_GIFT_Pollination, by.x = c('new_source_taxonID'), by.y= c('GBIF_taxonID'),all.x=T,allow.cartesian=TRUE)
harmonised_comp_rearr$glonaf_flag[harmonised_comp_rearr$new_source_taxonID %in% GBIF_GloNAF_uniq] <- 1
harmonised_comp_rearr$gift_flag[harmonised_comp_rearr$new_source_taxonID %in% GBIF_GIFT_Pollination$GBIF_taxonID] <- 1

full_df <- mutate(harmonised_comp_rearr,
                  infoFlags = case_when(
                    glonaf_flag == 1 & gift_flag == 1 ~ "GIFT & GloNAF",
                    glonaf_flag == 1 & is.na(gift_flag) ~ "GloNAF only",
                    is.na(glonaf_flag) & gift_flag == 1 ~ "GIFT only",
                    is.na(glonaf_flag) & is.na(gift_flag) ~ "NA")
)

#glonaf_sub <- full_df[work_ID!="NA" & glonaf_flag==1]
#glonaf_sub2 <- full_df[!is.na(infoFlags)]
#remove records with no species name or missing location
glonaf_sub2 <- subset(full_df,infoFlags!="NA" & new_source_specificEpithet != "" & decimalLongitude != "NA" & decimalLatitude!="NA")

top100_int <- fread("Top100_BioticInteractions_withGloNAF_GIFT_family_GIDEAS.csv")
top100_EICAT <- fread("Top100_Impacts_withGloNAF_GIFT_family.csv")

sf::sf_use_s2(FALSE)
test_species <- "Ailanthus altissima"
test_species <- ""


datalist = list()

for (test_species in top100_EICAT$Species) {

test_species_tab <-linked_tbl3[linked_tbl3$GBIF_species==test_species,]
test_species_int <-glonaf_sub2[glonaf_sub2$new_sourceName==test_species,]



# 4. Convert the data frame into a spatial sf object
# WGS84 (EPSG: 4326) is the standard coordinate system for GPS lat/long
pts_sf <- st_as_sf(test_species_int, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

# 5. Match the Coordinate Reference Systems (CRS)
# Both datasets must use the exact same CRS to intersect properly
test_species_tab_shape <- st_as_sf(test_species_tab[,c("OBJIDsic","IDregion","ISOcountry","LAT","LON","name.y","status","SRichness","Area_SqKm","geometry")])

# 6. Perform the spatial intersection (Spatial Join)
# This joins the shapefile attributes to the points they intersect
intersected_pts <- st_join(pts_sf, test_species_tab_shape, join = st_intersects)

#ggplot() + geom_sf(data=my_sf,fill="grey") + geom_sf(data=test_species_tab_shape,aes(fill = status), color = "white") + geom_sf(data=intersected_pts, size=1,alpha=0.5) + theme_void() + labs(title=test_species)
#ggsave(sprintf("./GloNAF_Ranges_x_GloBI_Impact_plots/%s.pdf", test_species), width = 12, height = 8)
datalist[[test_species]] <- intersected_pts
}

big_data = do.call(rbind, datalist)

#eventDate

intersected_pts_summary <- big_data %>% group_by(new_sourceName , status) %>% summarise(n = n())

ggplot(data=intersected_pts_summary, aes(y=new_sourceName, x=n, fill = status)) +labs(x="GloBI Interaction Frequency",y="Species") + geom_col() + xlim(0,5000)




sub <- subset(full_df, decimalLongitude == "NA" & decimalLatitude=="NA")
localities <- as.data.frame(table(sub$localityName))
localities %>% arrange(desc(Freq)) %>% slice(1:10)