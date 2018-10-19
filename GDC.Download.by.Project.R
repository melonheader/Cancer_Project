######Download GDC data by Project

#Attach libraries
require(TCGAbiolinks, quietly = TRUE)
require(tidyverse, quietly = TRUE)
require(SummarizedExperiment, quietly = TRUE)

#Select projects to download
Projects <- readline(prompt = "Projects to download: ")
print(paste("Selected projects: ", Projects, sep = ""))
Projects <- strsplit(Projects, ", ", fixed = TRUE) %>% 
  unlist() %>% 
  as.vector()

#Check directory
if(!dir.exists("GDCprocessed")){
  dir.create("GDCprocessed")
  print("Output directory GDCprocessed - created")
} else {
  print("Output directory GDCprocessed already exists")
}

#Write Data (HTSeq - Counts) for separate projects for later use
for(i in 1:length(Projects)) {
  #Prepare query
  Query <- GDCquery(project = Projects[i], 
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    workflow.type = "HTSeq - Counts",
                    sample.type = c("Primary solid Tumor","Solid Tissue Normal")
  )
  #Download data (if required)
  GDCdownload(Query)
  
  #Prepare matrix
  Exp.Query <- GDCprepare(query = Query, save = FALSE)
  Matr.Query <- assay(Exp.Query, "HTSeq - Counts")
  GeneID <- data.frame(GeneID = rownames(Matr.Query))
  
  #Write table
  write.table(Matr.Query, paste(Projects[i], ".Counts.tsv", sep = ""), 
              col.names = TRUE, 
              row.names = FALSE, 
              sep = "\t", 
              quote = FALSE)
  cat(paste(Projects[i], ".Counts.tsv is written to GDCprocessed"), sep = "")
} 

#Write universal geneIDS
write.table(GeneID, paste("TCGA", ".GeneIDs.tsv", sep = ""), 
            col.names = TRUE, 
            row.names = FALSE, 
            sep = "\t", 
            quote = FALSE)

#Clean environment
rm(Exp.Query, GeneID, Matr.Query, Query, i, Projects)
