import csv
import gzip
#from tqdm import tqdm
from rapidfuzz import process
import re
import Levenshtein
import sys

def process_taxon(sourceTaxonName):
    sourceTaxonName = sourceTaxonName.translate(str.maketrans({' ': ' ', '.': ' '}))
    sourceTaxonName_array = sourceTaxonName.split()
    if (len(sourceTaxonName_array) == 1):
        sName = sourceTaxonName_array[0]
    else:
        sName = sourceTaxonName_array[0] + ' ' + sourceTaxonName_array[1]

    if (sName in taxa_dict):
        # print('found in taxa_dict',taxa_dict[sName])
        taxonID, datasetID, parentNameUsageID, acceptedNameUsageID, originalNameUsageID, scientificName, scientificNameAuthorship, canonicalName, genericName, specificEpithet, infraspecificEpithet, taxonRank, nameAccordingTo, namePublishedIn, taxonomicStatus, nomenclaturalStatus, taxonRemarks, kingdom, phylum, class_clade, order, family, genus = \
        taxa_dict[sName]
        if (acceptedNameUsageID != ''):
            # print("synonym substitution performed before ",scientificName,":",canonicalName)
            if (acceptedNameUsageID in accepted_dict):
                taxonID, datasetID, parentNameUsageID, acceptedNameUsageID, originalNameUsageID, scientificName, scientificNameAuthorship, canonicalName, genericName, specificEpithet, infraspecificEpithet, taxonRank, nameAccordingTo, namePublishedIn, taxonomicStatus, nomenclaturalStatus, taxonRemarks, kingdom, phylum, class_clade, order, family, genus = \
                accepted_dict[acceptedNameUsageID]
                # print("synonym substitution performed")
                # print("synonym substitution performed after", scientificName,":", canonicalName)
            # else:
            # print(sName,acceptedNameUsageID)

    elif (sourceTaxonName_array[0] in taxa_dict):
        # print('found genus in taxa_dict', taxa_dict[sourceTaxonName_array[0]])
        taxonID, datasetID, parentNameUsageID, acceptedNameUsageID, originalNameUsageID, scientificName, scientificNameAuthorship, canonicalName, genericName, specificEpithet, infraspecificEpithet, taxonRank, nameAccordingTo, namePublishedIn, taxonomicStatus, nomenclaturalStatus, taxonRemarks, kingdom, phylum, class_clade, order, family, genus = \
            taxa_dict[sourceTaxonName_array[0]]
        if (acceptedNameUsageID != '' and acceptedNameUsageID in accepted_dict):
            # print("synonym genus substitution performed before ", scientificName, canonicalName)
            taxonID, datasetID, parentNameUsageID, acceptedNameUsageID, originalNameUsageID, scientificName, scientificNameAuthorship, canonicalName, genericName, specificEpithet, infraspecificEpithet, taxonRank, nameAccordingTo, namePublishedIn, taxonomicStatus, nomenclaturalStatus, taxonRemarks, kingdom, phylum, class_clade, order, family, genus = \
                accepted_dict[acceptedNameUsageID]
            # print("synonym substitution performed")
            # print("synonym genus substitution performed after", scientificName, canonicalName)

    else:
        taxonID, genericName, specificEpithet, kingdom, phylum, class_clade, order, family = "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"

        matches = process.extract(sName, taxa, limit=1)
        dist = Levenshtein.distance(matches[0][0], sName)
        # Protoeces hawaiiensis [('Proctoeces hawaiiensis', 97.67441860465115, 1051507)] 21
        if (dist < 2):
            #print(matches)
            #print(dist, 'close', sName, matches[0][0])
            taxonID, datasetID, parentNameUsageID, acceptedNameUsageID, originalNameUsageID, scientificName, scientificNameAuthorship, canonicalName,\
            genericName, specificEpithet, infraspecificEpithet, taxonRank,\
            nameAccordingTo, namePublishedIn, taxonomicStatus, nomenclaturalStatus,\
            taxonRemarks, kingdom, phylum, class_clade, order, family, genus = taxa_dict[matches[0][0]]

        else:
            #print('not', sName)
            taxonID, genericName, specificEpithet, kingdom, phylum, class_clade, order, family = "NA","NA","NA","NA","NA","NA","NA","NA"
    return([taxonID, genericName, specificEpithet, kingdom, phylum, class_clade, order, family])


def sliceindex(x):
    i = 0
    for c in x:
        if c.isalpha():
            i = i + 1
            return i
        i = i + 1

def upperfirst(x):
    i = sliceindex(x)
    return x[:i].upper() + x[i:]

GBIF = sys.argv[1]
taxa_dict= dict()
taxa = list()
accepted_dict = dict()
#non_accepted_dict = dict()
with open(GBIF) as f:
    header = f.readline()
    header2 = header.rstrip("\n")
    header2_sep = header2.split("\t")
    #result = ['sourceGBIF_' + element for element in header2_sep]
    #print(result)
    for line in f:
        line2 = line.rstrip("\n")
        sep = line2.split("\t")

        taxonID, datasetID, parentNameUsageID, acceptedNameUsageID, originalNameUsageID, scientificName, scientificNameAuthorship,\
        canonicalName, genericName, specificEpithet, infraspecificEpithet, taxonRank, nameAccordingTo, namePublishedIn, \
        taxonomicStatus, nomenclaturalStatus, taxonRemarks, kingdom, phylum, class_clade, order, family,	genus = sep
        #print(scientificName, kingdom, phylum)


        ###process strings

        scientificName = scientificName.replace('.', ' ')
        canonicalName = canonicalName.replace('.', ' ')
        scientificName = upperfirst(scientificName)
        canonicalName = upperfirst(canonicalName)
        #print(scientificName)
        #print(canonicalName)
        if(taxonomicStatus=="accepted"):
            accepted_dict[taxonID] = sep

        if(canonicalName==""):
            sciName_array = scientificName.split(" ")
            if(len(sciName_array)==1):
                taxa_dict[sciName_array[0]] = sep
                taxa.append(sciName_array[0])
            else:
                taxa_dict[sciName_array[0]+" "+sciName_array[1]] = sep
                taxa.append(sciName_array[0]+" "+sciName_array[1])
                #print(sciName_array[0]+" "+sciName_array[1])
        else:
            canName_array = canonicalName.split(" ")
            if (len(canName_array) == 1):
                taxa_dict[canName_array[0]] = sep
                taxa.append(canName_array[0])
            else:
                taxa_dict[canName_array[0] + " " + canName_array[1]] = sep
                taxa.append(canName_array[0] + " " + canName_array[1])
                #print(canName_array[0] + " " + canName_array[1])


# sourceTaxonId,sourceTaxonIds,sourceTaxonName,sourceTaxonRank,sourceTaxonPathNames,sourceTaxonPathIds,
# sourceTaxonPathRankNames,sourceTaxonSpeciesName,sourceTaxonSpeciesId,sourceTaxonSubgenusName,
# sourceTaxonSubgenusId,sourceTaxonGenusName,sourceTaxonGenusId,sourceTaxonFamilyName,sourceTaxonFamilyId,
# sourceTaxonOrderName,sourceTaxonOrderId,sourceTaxonClassName,sourceTaxonClassId,sourceTaxonPhylumName,
# sourceTaxonPhylumId,sourceTaxonKingdomName,sourceTaxonKingdomId,sourceId,sourceOccurrenceId,
# sourceInstitutionCode,sourceCollectionCode,sourceCatalogNumber,sourceBasisOfRecordId,sourceBasisOfRecordName,
# sourceLifeStageId,sourceLifeStageName,sourceBodyPartId,sourceBodyPartName,sourcePhysiologicalStateId,
# sourcePhysiologicalStateName,sourceSexId,sourceSexName,interactionTypeName,interactionTypeId,targetTaxonId,
# targetTaxonIds,targetTaxonName,targetTaxonRank,targetTaxonPathNames,targetTaxonPathIds,targetTaxonPathRankNames,
# targetTaxonSpeciesName,targetTaxonSpeciesId,targetTaxonSubgenusName,targetTaxonSubgenusId,targetTaxonGenusName,
# targetTaxonGenusId,targetTaxonFamilyName,targetTaxonFamilyId,targetTaxonOrderName,targetTaxonOrderId,
# targetTaxonClassName,targetTaxonClassId,targetTaxonPhylumName,targetTaxonPhylumId,targetTaxonKingdomName,
# targetTaxonKingdomId,targetId,targetOccurrenceId,targetInstitutionCode,targetCollectionCode,targetCatalogNumber,
# targetBasisOfRecordId,targetBasisOfRecordName,targetLifeStageId,targetLifeStageName,targetBodyPartId,targetBodyPartName,
# targetPhysiologicalStateId,targetPhysiologicalStateName,targetSexId,targetSexName,decimalLatitude,decimalLongitude,
# localityId,localityName,eventDate,argumentTypeId,referenceCitation,referenceDoi,referenceUrl,sourceCitation,
# sourceNamespace,sourceArchiveURI,sourceDOI,sourceLastSeenAtUnixEpoch

column_name = sys.argv[3]

#GloBI = "/Users/Sam/Desktop/PostDoc_Ecology/interactions_plantae.csv.gz"
GloBI = sys.argv[2]
with gzip.open(GloBI, mode='rt') as G:
#with open(GloBI) as G:
    header_G = G.readline()
    #globi_header = csv.reader(header_G, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True)
    #print(globi_header)
    globi_header2 = header_G.rstrip("\n")
    globi_header2_sep = globi_header2.split(",")
    #print(globi_header2_sep)

    source_header = ["\"GBIF_source_taxonID\"", "\"GBIF_source_genericName\"", "\"GBIF_source_specificEpithet\"", "\"GBIF_source_kingdom\"", "\"GBIF_source_phylum\"", "\"GBIF_source_class\"", "\"GBIF_source_order\"", "\"GBIF_source_family\""]

    header_list = source_header + globi_header2_sep
    print(','.join(header_list))
    for l in csv.reader(G, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True):
        #print(len(l))
        res = dict(zip(globi_header2_sep, l))
        #work_ID,work_species,work_author,trait_value_361,trait_value_362,trait_value_363,agreement_361,agreement_362,agreement_363,n_361,n_362,n_363,references_361,references_362,references_363 = l

        #sourceTaxonName_array = sourceTaxonName.split(" ")
        gbif_source = process_taxon(res[column_name])
        output_list = gbif_source + l

        print(', '.join(f'"{w}"' for w in output_list))
#

# with open(GloBI) as f:
#     f.readline()
#     for line in tqdm(f):
#         line2 = line.rstrip("\n")
#         sep = line2.split(",")
#         sourceTaxonId, sourceTaxonIds, sourceTaxonName, sourceTaxonRank, sourceTaxonPathNames, sourceTaxonPathIds,\
#         sourceTaxonPathRankNames,sourceTaxonSpeciesName,sourceTaxonSpeciesId,sourceTaxonSubgenusName,\
#         sourceTaxonSubgenusId,sourceTaxonGenusName,sourceTaxonGenusId,sourceTaxonFamilyName,sourceTaxonFamilyId,\
#         sourceTaxonOrderName,sourceTaxonOrderId,sourceTaxonClassName,sourceTaxonClassId,sourceTaxonPhylumName,\
#         sourceTaxonPhylumId,sourceTaxonKingdomName,sourceTaxonKingdomId,sourceId,sourceOccurrenceId,\
#         sourceInstitutionCode,sourceCollectionCode,sourceCatalogNumber,sourceBasisOfRecordId,sourceBasisOfRecordName,\
#         sourceLifeStageId,sourceLifeStageName,sourceBodyPartId,sourceBodyPartName,sourcePhysiologicalStateId,\
#         sourcePhysiologicalStateName,sourceSexId,sourceSexName,interactionTypeName,interactionTypeId,targetTaxonId,\
#         targetTaxonIds,targetTaxonName,targetTaxonRank,targetTaxonPathNames,targetTaxonPathIds,targetTaxonPathRankNames,\
#         targetTaxonSpeciesName,targetTaxonSpeciesId,targetTaxonSubgenusName,targetTaxonSubgenusId,targetTaxonGenusName,\
#         targetTaxonGenusId,targetTaxonFamilyName,targetTaxonFamilyId,targetTaxonOrderName,targetTaxonOrderId,\
#         targetTaxonClassName,targetTaxonClassId,targetTaxonPhylumName,targetTaxonPhylumId,targetTaxonKingdomName,\
#         targetTaxonKingdomId,targetId,targetOccurrenceId,targetInstitutionCode,targetCollectionCode,targetCatalogNumber,\
#         targetBasisOfRecordId,targetBasisOfRecordName,targetLifeStageId,targetLifeStageName,targetBodyPartId,targetBodyPartName,\
#         targetPhysiologicalStateId,targetPhysiologicalStateName,targetSexId,targetSexName,decimalLatitude,decimalLongitude,\
#         localityId,localityName,eventDate,argumentTypeId,referenceCitation,referenceDoi,referenceUrl,sourceCitation,\
#         sourceNamespace,sourceArchiveURI,sourceDOI,sourceLastSeenAtUnixEpoch = sep
#         if(sourceTaxonName in taxa_dict):
#             print('found')
#         else:
#             print('not')
#
#

