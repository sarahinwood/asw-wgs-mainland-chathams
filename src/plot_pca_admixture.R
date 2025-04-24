library(tidyverse)
library(data.table)
library(viridis)
library(ggplot2)

pca_pruned <- fread('output/04_plink/ld_pruned/ld_pruned_snps_plink_pca.eigenvec')
eigenval_pruned <- scan('output/04_plink/ld_pruned/ld_pruned_snps_plink_pca.eigenval')

# set names
names(pca_pruned)[1] <- "sample_name"

pca_pruned$Location <- tstrsplit(pca_pruned$sample_name, "_", keep=1)
pca_pruned$Location <- ifelse(pca_pruned$Location=="Chatham", paste("Chatham Island"), 
                              ifelse(pca_pruned$Location=="Stewart", paste("Stewart Island"), paste(pca_pruned$Location)))

pca_pruned$Location <- factor(pca_pruned$Location, levels=c("Chatham Island", "Stewart Island", "Fortrose", "Lincoln"))

# first convert to percentage variance explained
pve_pruned <- data.frame(PC = 1:20, pve = eigenval_pruned/sum(eigenval_pruned)*100)

# calculate the cumulative sum of the percentage variance explained
cumsum(pve_pruned$pve)

# plot pca_pruned - Location
ggplot(pca_pruned, aes(PC1, PC2, colour=Location))+ 
  geom_point(size = 3, alpha=0.7)+
  scale_colour_viridis(discrete=T)+ 
  coord_equal()+
  theme_light()+
  xlab(paste0("PC1 (", signif(pve_pruned$pve[1], 3), "%)"))+
  ylab(paste0("PC2 (", signif(pve_pruned$pve[2], 3), "%)"))

pca_pruned$sample_name <- factor(pca_pruned$sample_name, levels=pca_pruned$sample_name[order(pca_pruned$Location, pca_pruned$sample_name, decreasing=F)])
long_pca_pruned <- melt(pca_pruned)

# Create a named vector for facet labels
pc_labels_pruned <- setNames(
  paste0("PC", pve_pruned$PC, " (", signif(pve_pruned$pve, 3), "%)"), 
  paste0("PC", pve_pruned$PC)
)

# plot all pcas_pruned - location
ggplot(long_pca_pruned, aes(sample_name, value, colour=Location))+
  geom_point(size = 3, alpha=0.7)+
  scale_colour_viridis(discrete=T)+ 
  theme_light()+
  facet_wrap(~variable, labeller = labeller(variable = pc_labels_pruned)) +
  labs(x = NULL,  # Remove x-axis label
       y = "PC value") +  # Change y-axis label
  theme(axis.text.x = element_blank(),  # Remove x-axis text (sample names)
        axis.ticks.x = element_blank())  # Remove x-axis ticks
