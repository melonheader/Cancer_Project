######Prepare matrix from downloaded GDC data

#Attach libraries
require(TCGAbiolinks)
require(tidyverse, quietly = TRUE)

#Select projects to prepare
Selection <- readline(prompt = "Select projects to prepare: ")
cat(paste("Selected projects: ", Selection, sep = ""))
Selection.M <- strsplit(Selection, ", ", fixed = TRUE) %>% 
  unlist()
Selection <- strsplit(Selection, ", ", fixed = TRUE) %>% 
  unlist() %>%
  paste(., ".Counts.tsv", sep = "")

#Read GDC data into one dataframe
Data.Path <- "GDCprocessed"
Data <- Selection %>%
  map_dfc(~read_tsv(file.path(Data.Path, .)))

#Attach gene id
GeneIDs <- read_tsv(file.path(Data.Path, dir(Data.Path, "*IDs.tsv")))
Data <- Data %>% 
  mutate(Gene_ID = GeneIDs$GeneID) %>%
  select(., Gene_ID, everything())
rm(GeneIDs)

#Attach non-GDC projects
Filepaths <- character(0)
repeat {
  Inquiry <- readline(prompt = "Attach non-GDC data? ")
  if(Inquiry == "y") {
    Data.Path.1 <- readline(prompt = "Specify path to data: ")
    if(!dir.exists(file.path(Data.Path.1))){
      print("Specified path is incorrect")
    } else {
      Filepaths <- c(Filepaths, Data.Path.1)
      print(paste(Data.Path.1, " added", sep = ""))
    }
  } else {
    print("Proceed.....")
    #Read selected non-GDC data
    Data.List <- lapply(Filepaths, read_tsv)
    names(Data.List) <- basename(Filepaths)
    list2env(Data.List, .GlobalEnv)
    rm(Data.List)
    print(paste(basename(Filepaths), "are added succesfully", sep = ""))
    break
  }
}


#Get metadata
print(paste("Collecting metadata for ", Selection.M, " Project", sep = ""))
Query <- GDCquery(project = Selection.M, 
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts",
                  sample.type = c("Primary solid Tumor","Solid Tissue Normal")
)
#Which samples are primary solid tumor
Data.Sm.TP <- TCGAquery_SampleTypes(getResults(Query, cols="cases"),"TP") %>% 
  data.frame("Sm.Barcode" = ., 
             "Sm.Type" = rep("Primary Tumor", length(.)))
#Which samples are normal solid tissue
Data.Sm.NT <- TCGAquery_SampleTypes(getResults(Query, cols="cases"),"NT") %>% 
  data.frame("Sm.Barcode" = ., 
             "Sm.Type" = rep("Normal Tissue", length(.)))
#Merge
Data.Sm <- rbind(Data.Sm.TP, Data.Sm.NT) %>%
  mutate(bcr_patient_barcode = substring(Sm.Barcode, 1, 12)) %>%
  select(bcr_patient_barcode, everything())
rm(Data.Sm.TP, Data.Sm.NT)

#Get clinical data
Data.Clin <- data.frame()
for (i in 1:length(Selection.M)) {
  Data.Clin.Sep <- GDCquery_clinic(project = Selection.M[i], "clinical")
  Data.Clin <- rbind(Data.Clin, 
                     Data.Clin.Sep[, c("bcr_patient_barcode", "disease", 
                                       "gender", "tumor_stage", "tumor_grade", 
                                       "age_at_diagnosis","vital_status", 
                                       "days_to_death")])
  rm(Data.Clin.Sep, i)
}

#Merge TCGA sample and clinical metadata
Data.Meta <- full_join(Data.Clin, Data.Sm, by = "bcr_patient_barcode") %>%
  select(bcr_patient_barcode, Sm.Barcode, Sm.Type, everything())
rm(Data.Clin, Data.Sm)
#Check directory
if(!dir.exists("DataClinical")){
  dir.create("DataClinical")
  print("Output directory DataClinical - created")
} else {
  print("Output directory DataClinical already exists")
}

#Write metadata for selected projects
write_tsv(Data.Meta, 
          file.path("DataClinical", paste("Data.Meta.Tcga.", Sys.Date(), ".tsv", sep = "")),
          col_names = TRUE)

print(paste("Metadata for ", Selection.M, 
            " was written to DataClinical/", 
            sep = ""))
rm(Selection.M, Data.Clin, Data.Sm)

#Create and standartize clinical data for non-GDC data
Data.Clin.gtex <- read_csv("Data/gtex.Annotation.csv", 
                             col_names = TRUE) %>%
  filter(., SMATSSCR <= 1) %>%
  transmute(bcr_patient_barcode = SAMPID,
            Sm.Barcode = SAMPID,
            Sm.Type = SMTS,
            disease = SMTSD, 
            gender = rep(NA, dim(.)[1]), 
            tumor_stage = rep("Normal Tissue", dim(.)[1]), 
            tumor_grade = rep("Normal Tissue", dim(.)[1]), 
            age_at_diagnosis = rep(NA, dim(.)[1]), 
            vital_status = rep(NA, dim(.)[1]), 
            days_to_death = rep(NA, dim(.)[1])
            )
Data.Clin.De <- tibble(bcr_patient_barcode = names(DE.Counts.tsv[2:dim(DE.Counts.tsv)[2]]),
                       Sm.Barcode = names(DE.Counts.tsv[2:dim(DE.Counts.tsv)[2]]),
                       Sm.Type = names(DE.Counts.tsv[2:dim(DE.Counts.tsv)[2]]),
                       disease = rep("Definitive Endoderm", dim(DE.Counts.tsv)[2] - 1), 
                       gender = rep(NA, dim(DE.Counts.tsv)[2] - 1), 
                       tumor_stage = rep("Definitive Endoderm", dim(DE.Counts.tsv)[2] - 1), 
                       tumor_grade = rep("Definitive Endoderm", dim(DE.Counts.tsv)[2] - 1), 
                       age_at_diagnosis = rep(NA, dim(DE.Counts.tsv)[2] - 1), 
                       vital_status = rep(NA, dim(DE.Counts.tsv)[2] - 1), 
                       days_to_death = rep(NA, dim(DE.Counts.tsv)[2] - 1)
                       )

#Write metadata for selected projects
write_tsv(Data.Clin.gtex, 
          file.path("DataClinical", paste("Data.Clin.gtex.", Sys.Date(), ".tsv", sep = "")),
          col_names = TRUE)
write_tsv(Data.Clin.De, 
          file.path("DataClinical", paste("Data.Clin.De.", Sys.Date(), ".tsv", sep = "")),
          col_names = TRUE)
print(paste("Metadata for ", basename(Filepaths), 
            " was written to DataClinical/", 
            sep = ""))
rm(Data.Clin.De, Data.Clin.gtex)

###Dalshe ruchkami
gtex.Counts.tsv$Description <- NULL
names(gtex.Counts.tsv)[1] <- "Gene_ID"
gtex.Counts.tsv$Gene_ID <- gsub("\\..*","",gtex.Counts.tsv$Gene_ID)

#Merge all data
Data.Merged <- inner_join(Data, DE.Counts.tsv, by = "Gene_ID") %>%
  inner_join(., gtex.Counts.tsv, by = "Gene_ID") %>%
  na.omit()
rm(DE.Counts.tsv, gtex.Counts.tsv, Data)

#Write Data
if(!dir.exists("DataTemp")){
  dir.create("DataTemp")
  print("Output directory DataTemp - created")
} else {
  print("Output directory DataTemp already exists")
}
write_tsv(Data.Merged, file.path("DataTemp", 
                                 paste("Data.Merged.", Sys.Date(), ".tsv", 
                                       sep = "")), 
          col_names = TRUE)
#Clean
rm(list=ls(all.names = TRUE))
.rs.restartR()
