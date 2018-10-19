######Normalize expression data


#Plot medians (dplyr)
#Median.log.Cpm <- Cpm.log %>% 
#  as.tibble() %>%
#  rowwise() %>% 
#  transmute(
#    Median = median(c(!!! rlang::syms(names(.))))
#  )
#ggplot(Median.log.Cpm, aes(x = Median)) +
#  geom_histogram() +
#  geom_vline(aes(xintercept = Cutoff), colour = "red")


#Attach libraries
require(TCGAbiolinks)
require(tidyverse, quietly = TRUE)
require(edgeR)
require(ggplot2)

#Get data
Data.Merged <- read_tsv("DataTemp/Data.Merged.2018-10-11.tsv", 
                        col_names = TRUE)
Matr.Comp <- Data.Merged %>% 
  select(-Gene_ID) %>%
  na.omit() %>%
  as.data.frame()
GeneIDs <- Data.Merged %>%
  select(Gene_ID)
rm(Data.Merged)

#Quality control (log cpm cutoff)
Cpm.Log <- cpm(Matr.Comp, log = TRUE)
Median.log2.cpm <- apply(Cpm.Log, 1, median)
hist(Median.log2.cpm)
Expr.Cutoff <- -4
abline(v = Expr.Cutoff, col = "red", lwd = 3)
sum(Median.log2.cpm > Expr.Cutoff)
rm(Cpm.Log)
rownames(Matr.Comp) <- GeneIDs$Gene_ID
counts <- Matr.Comp[Median.log2.cpm > Expr.Cutoff, ]
rm(Matr.Comp)

#Normalize with edgeR
y <- DGEList(counts, genes = rownames(counts))
y <- calcNormFactors(y, method = "TMM")
Normalized.CPM <- cpm(y, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 0.1)
Normalized.CPM <- Normalized.CPM %>%
  as.tibble(.) %>%
  mutate(GeneID = y$genes$genes) %>%
  select(GeneID, everything())

#Checkpoint
write_tsv(Normalized.CPM, 
            file.path("DataTemp", 
                      paste("Normalized.Counts.", Sys.Date(), ".tsv", sep = "")),
          append = FALSE)
#Clean
rm(list=ls(all.names = TRUE))
.rs.restartR()