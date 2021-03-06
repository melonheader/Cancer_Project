######Filter and add clinical data

#Attach libraries
#require(TCGAbiolinks)
require(tidyverse)
require(ggplot2)
require(factoextra)
require(FactoMineR)

##PCA
#Check components
Data <- Data[!is.na(Data$Sm.Type), ]
#Comp.Exp.Ann[, -c(1:10)] <- lapply(Comp.Exp.Ann[, -c(1:10)], function(x) as.numeric(as.character(x)))
Data.Pca <- PCA(Data[, -c(1:10)],  graph = FALSE, scale.unit = FALSE)
fviz_contrib(Data.Pca, choice = "var", axes = 1, top = 30)
fviz_eig(Data.Pca, addlabels = TRUE, ylim = c(0, 50))

plot_var <- fviz_pca_var(Data.Pca, col.var = "contrib",
                         ggtheme = theme_minimal(), gradient.cols = c("blue", "red")
)


#Tumor stage PCA

#Days to death PCA with gradient code
require(ggplot2)
Res.Pca <- prcomp(Data[, -c(1:10)],
                  center = TRUE,
                  scale. = FALSE
)

#Annotate to plot gradient coloration
Results.Pca <- as.data.frame(Res.Pca$x)
Results.Pca <- Results.Pca %>%
  as.tibble() %>%
  mutate(sample_type = Data$Sm.Type, 
         tumor_type = Data$disease, 
         tumor_stage = Data$tumor_stage, 
         days_to_death = Data$days_to_death) %>%
  select(sample_type, tumor_type, tumor_stage, days_to_death, everything() )

#Separate data with days_to_death
Results.Pca.surv <- Results.Pca[!is.na(Results.Pca$days_to_death), ]
#Cut below 2000 days_to_death
Results.Pca.surv <- Results.Pca.surv[Results.Pca.surv$days_to_death < 2100, ]
Results.Pca.nsurv <- Results.Pca[is.na(Results.Pca$days_to_death), ]


ggplot() +
  stat_ellipse(data = Results.Pca,
               aes(x = Results.Pca$PC1, y = Results.Pca$PC2, 
                   fill = Results.Pca$sample_type),
               type = "norm",
               geom = "polygon",
               level = 0.8,
               alpha = 0.05
  ) +
  geom_point(data = Results.Pca.surv,
             aes(x = Results.Pca.surv$PC1, y = Results.Pca.surv$PC2, 
                 col = Results.Pca.surv$days_to_death), 
             alpha = 0.8
  ) +
  geom_point(data = Results.Pca.nsurv,
             aes(x = Results.Pca.nsurv$PC1, y = Results.Pca.nsurv$PC2), 
             alpha = 0.05
  ) +
scale_colour_gradient2(low = "darkslateblue",mid = "white", high = "brown3"
) + 
scale_fill_manual(values=c("chocolate2", "springgreen4", "red3")
) +
labs(
  x = "PC1 (24%)",
  y = "PC2 (13.5%)",
  colour = "Days to Death",
  title = "PCA of various tumors (LUAD, LUSC, LIHC, GBB, LGG, SARC)",
  fill = "Sample Type"
)

##t-SNE
## Executing tsne
require(Rtsne)
Tsne <- Rtsne(Data[, -c(1:10)], 
              dims = 2, 
              perplexity = 85, 
              verbose = TRUE, 
              max_iter = 1000
)
#exeTimeTsne<- system.time(Rtsne(train[,-1], dims = 2, perplexity=30, verbose=TRUE, max_iter = 500))
Tsne.Plot <- data.frame(x = Tsne$Y[ ,1], y = Tsne$Y[ ,2], 
                        tumor_stage = Data$tumor_stage, 
                        disease = Data$disease,
                        days_to_death = Data$days_to_death,
                        Sample_Type = Data$Sm.Type
)
ggplot(Tsne.Plot) +
  stat_ellipse(aes(x = x, y = y, 
                   fill = Tsne.Plot$Sample_Type),
               type = "norm",
               geom = "polygon",
               level = 0.8,
               alpha = 0.2
  ) +
  geom_point(aes(x = x, y = y,
                 col = disease,
                 alpha = 0.7
  )
  ) +
  #scale_colour_manual(values = c("turquoise4", "red3")
  #) + 
  #scale_fill_manual(values = c("chocolate2", "springgreen4")
  #) +
  labs(
    x = "Dim 1",
    y = "Dim 2",
    colour = "Tumor type",
    title = "t-SNE LUAD, LUSC",
    fill = "Sample Type"
  )


#Save plot button
#Save gradient
if(!dir.exists("Plots")){
  dir.create("Plots")
  print("Output directory Plots - created")
} else {
  print("Output directory Plots already exists")
}

ggsave(file.path("Plots.", paste("PCA_gradient.", Sys.Date(), ".png", sep = "")), 
       plot = plot_sample, 
       device = "png", 
       units = "mm",
       width = 300,
       heigh = 160,
       dpi = 500)

#Save PCA
#Save t-sne
if(!dir.exists("Plots")){
  dir.create("Plots")
  print("Output directory Plots - created")
} else {
  print("Output directory Plots already exists")
}
ggsave("Plots/tSNE_LUAD_LUSC.png", 
       plot = last_plot(), 
       device = "png", 
       units = "mm",
       width = 300,
       heigh = 160,
       dpi = 500)

