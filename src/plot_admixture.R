library(tidyverse)
library(ggplot2)
library(forcats)
library(data.table)
library(viridis)

## which K value best?
cv_res <- fread('output/04_plink/ld_pruned/admixture/admixture_cv_res.out', header=F)
setorder(cv_res, V4)
cv_res$K_numeric <- as.numeric(gsub("\\D", "", cv_res$V3))

# Plot
ggplot(cv_res, aes(x = K_numeric, y = V4)) +
  geom_line() +
  geom_point() +
  labs(x = "Number of Clusters (K)", y = "Cross-Validation Error") +
  theme_minimal()

## k=2 ##

admixture_2 <- read_delim("output/04_plink/ld_pruned/admixture/admixture.2.Q",
                          col_names = paste0("Q",seq(1:2)),
                          delim=" ")
admixture_2$k <- "2"

sample_names_table <- fread("output/04_plink/ld_pruned/ld_pruned_snps_plink_pca.fam", header=F)
sample_names <- sample_names_table[,2]
admixture_2$sample <- sample_names$V2

long_admix2 <- gather(admixture_2, Q, value, -sample,-k)
## get location names
long_admix2 <- data.table(long_admix2)
long_admix2$Location <- ifelse(grepl("Chatham", long_admix2$sample), paste("Chatham Island"),
                               ifelse(grepl("Stewart", long_admix2$sample), paste("Stewart Island"),
                                      ifelse(grepl("Lincoln", long_admix2$sample), paste("Lincoln"),
                                             ifelse(grepl("Fortrose", long_admix2$sample), paste("Fortrose"), paste(long_admix2$Location)))))
long_admix2$Location <- factor(long_admix2$Location, levels=c("Chatham Island", "Stewart Island", "Fortrose", "Lincoln"))

ggplot(long_admix2, aes(x=sample, y=value, fill=factor(Q)))+ 
  geom_col(width=1)+
  facet_grid(~ Location, scales = "free_x", space = "free_x", switch = "x") +  # Separate locations with gaps
  xlab("Population")+
  ylab("Ancestry")+
  scale_y_continuous(limits = c(0, 1))+
  theme_minimal()+
  scale_fill_viridis(discrete=T, option="C", direction=-1)+
  theme(axis.ticks.y = element_line(color = "black"),
        panel.grid.major = element_blank(),  # Remove major horizontal grid lines
        panel.grid.minor = element_blank(),  # Remove minor horizontal grid lines
        panel.spacing.x = unit(1, "mm"),  # Add space between locations
        axis.text.x = element_blank(),  # Hide sample names
        strip.placement = "outside",  # Moves the location labels below
        strip.text = element_text(size = 10, face = "bold"),
        strip.switch.pad.grid = unit(-6, "mm"),  # **Reduces space between graph and location labels**
        plot.margin = margin(5, 5, 5, 5), # **Reduce bottom margin**
        legend.position = "none")


## k=3 ##

admixture_3 <- read_delim("output/04_plink/ld_pruned/admixture/admixture.3.Q",
           col_names = paste0("Q",seq(1:3)),
           delim=" ")
admixture_3$k <- "3"

admixture_3$sample <- sample_names$V2


long_admix3 <- gather(admixture_3, Q, value, -sample,-k)
## get location names
long_admix3 <- data.table(long_admix3)
long_admix3$Location <- ifelse(grepl("Chatham", long_admix3$sample), paste("Chatham Island"),
                                     ifelse(grepl("Stewart", long_admix3$sample), paste("Stewart Island"),
                                                  ifelse(grepl("Lincoln", long_admix3$sample), paste("Lincoln"),
                                                               ifelse(grepl("Fortrose", long_admix3$sample), paste("Fortrose"), paste(long_admix3$Location)))))
long_admix3$Location <- factor(long_admix3$Location, levels=c("Chatham Island", "Stewart Island", "Fortrose", "Lincoln"))

ggplot(long_admix3, aes(x=sample, y=value, fill=factor(Q)))+ 
  geom_col(width=1)+
  facet_grid(~ Location, scales = "free_x", space = "free_x", switch = "x") +  # Separate locations with gaps
  xlab("Population")+
  ylab("Ancestry")+
  scale_y_continuous(limits = c(0, 1.0001))+
  theme_minimal()+
  scale_fill_viridis(discrete=T, option="C", direction=-1)+
  theme(axis.ticks.y = element_line(color = "black"),
        panel.grid.major = element_blank(),  # Remove major horizontal grid lines
        panel.grid.minor = element_blank(),  # Remove minor horizontal grid lines
        panel.spacing.x = unit(1, "mm"),  # Add space between locations
        axis.text.x = element_blank(),  # Hide sample names
        strip.placement = "outside",  # Moves the location labels below
        strip.text = element_text(size = 10, face = "bold"),
        strip.switch.pad.grid = unit(-6, "mm"),  # **Reduces space between graph and location labels**
        plot.margin = margin(5, 5, 5, 5), # **Reduce bottom margin**
        legend.position = "none")

## k=4 ##

admixture_4 <- read_delim("output/04_plink/ld_pruned/admixture/admixture.4.Q",
                          col_names = paste0("Q",seq(1:4)),
                          delim=" ")
admixture_4$k <- "4"

admixture_4$sample <- sample_names$V2


long_admix4 <- gather(admixture_4, Q, value, -sample,-k)
## get location names
long_admix4 <- data.table(long_admix4)
long_admix4$Location <- ifelse(grepl("Chatham", long_admix4$sample), paste("Chatham Island"),
                               ifelse(grepl("Stewart", long_admix4$sample), paste("Stewart Island"),
                                      ifelse(grepl("Lincoln", long_admix4$sample), paste("Lincoln"),
                                             ifelse(grepl("Fortrose", long_admix4$sample), paste("Fortrose"), paste(long_admix4$Location)))))
long_admix4$Location <- factor(long_admix4$Location, levels=c("Chatham Island", "Stewart Island", "Fortrose", "Lincoln"))

ggplot(long_admix4, aes(x=sample, y=value, fill=factor(Q)))+ 
  geom_col(width=1)+
  facet_grid(~ Location, scales = "free_x", space = "free_x", switch = "x") +  # Separate locations with gaps
  xlab("Population")+
  ylab("Ancestry")+
  scale_y_continuous(limits = c(0, 1.000001))+
  theme_minimal()+
  scale_fill_viridis(discrete=T, option="C", direction=1)+
  theme(axis.ticks.y = element_line(color = "black"),
        panel.grid.major = element_blank(),  # Remove major horizontal grid lines
        panel.grid.minor = element_blank(),  # Remove minor horizontal grid lines
        panel.spacing.x = unit(1, "mm"),  # Add space between locations
        axis.text.x = element_blank(),  # Hide sample names
        strip.placement = "outside",  # Moves the location labels below
        strip.text = element_text(size = 10, face = "bold"),
        strip.switch.pad.grid = unit(-6, "mm"),  # **Reduces space between graph and location labels**
        plot.margin = margin(5, 5, 5, 5), # **Reduce bottom margin**
        legend.position = "none")
