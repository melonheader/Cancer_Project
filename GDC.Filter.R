######Filter and add clinical data

#Attach libraries
#require(TCGAbiolinks)
require(tidyverse)
require(ggplot2)

#Read Normalized data
Data <- read_tsv(file.path("DataTemp", 
                           paste("Normalized.Counts.", Sys.Date(), ".tsv", sep = "")), 
                 col_names = TRUE)

#Load gene.set
Gene.Set.1 <- read_table("Markers/ESC_DP.csv", col_names = TRUE)

#Change EnsID to hgncID
require(biomaRt)
Mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
Gene.List <- getBM(filters = "ensembl_gene_id", 
                   attributes = c("ensembl_gene_id","hgnc_symbol"),
                   values = Data$GeneID, 
                   mart = Mart
)
rm(Mart)
detach("package:biomaRt", unload=TRUE)

Gene.List <- Gene.List %>%
  as.tibble() %>%
  filter(., hgnc_symbol != "")

#Filter gene set
Data <- Data %>%
  right_join(., Gene.List, by = c("GeneID" = "ensembl_gene_id")) %>%
  rename(., Gene_ID = hgnc_symbol) %>%
  select(-GeneID) %>%
  filter(Gene_ID %in% Gene.Set.1$GeneID) %>%
  select(Gene_ID, everything())

#Transpose
Data <- Data %>%
  gather(SampleID, Counts, 2:ncol(Data)) %>%
  spread(Gene_ID, Counts)

#Load metdata
Filepaths <- list.files("DataClinical")
Data.Clin.Path <- "DataClinical"
Data.Meta <- Filepaths %>%
  map_dfr(~read_tsv(file.path(Data.Clin.Path, .), 
                    col_types = cols(
                      bcr_patient_barcode = col_character(),
                      Sm.Barcode = col_character(),
                      Sm.Type = col_character(),
                      disease = col_character(),
                      gender = col_character(),
                      tumor_stage = col_character(),
                      tumor_grade = col_character(),
                      age_at_diagnosis = col_integer(),
                      vital_status = col_character(),
                      days_to_death = col_double()
                    )))


#Attach metadata
Data <- Data.Meta %>%
  right_join(., Data, by = c("Sm.Barcode" = "SampleID"))

#Standartize tumor stages
Data[Data$tumor_stage %in% c("stage i", "stage ia", "stage ib"), ]$tumor_stage <- 
  rep("Early stage", 
      length(Data[Data$tumor_stage %in% c("stage i", "stage ia", "stage ib"), ]$tumor_stage)
  )
Data[Data$tumor_stage %in% c("stage ii", "stage iia", "stage iib"), ]$tumor_stage <- 
  rep("Intermediate stage", 
      length(Data[Data$tumor_stage %in% c("stage ii", "stage iia", "stage iib"), ]$tumor_stage)
  )
Data[Data$tumor_stage %in% c("stage iii", "stage iiia", "stage iiib", "stage iiic", "stage iv", "stage iva", "stage ivb"), ]$tumor_stage <- 
  rep("Late stage", 
      length(Data[Data$tumor_stage %in% c("stage iii", "stage iiia", "stage iiib", "stage iiic", "stage iv", "stage iva", "stage ivb"), ]$tumor_stage)
  )

