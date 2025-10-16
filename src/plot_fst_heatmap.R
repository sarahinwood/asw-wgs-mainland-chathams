library(data.table)
library(tidyverse)
library(reshape2)
library(pheatmap)
library(viridis)
library(reshape2)

####################
## location based ##
####################

# Read the PLINK2 fst matrix
fst <- fread("output/04_plink/fst/Location/population_fst_results.fst.summary")

sample_table <- fread("data/sample_table.csv")
sample_to_location <- sample_table[,c(1,7)]

# Get unique population names
populations <- unique(c(fst$`#POP1`, fst$POP2))
# Initialize empty matrix
fst_matrix <- matrix(NA, nrow=length(populations), ncol=length(populations),
                     dimnames=list(populations, populations))

# Fill matrix using data from lower-triangle format
for (i in 1:nrow(fst)) {
  pop1 <- fst$`#POP1`[i]
  pop2 <- fst$POP2[i]
  fst_value <- fst$HUDSON_FST[i]
  
  fst_matrix[pop1, pop2] <- fst_value
  fst_matrix[pop2, pop1] <- fst_value  # Make symmetric
}

fst_matrix_ord <- fst_matrix[c(1,4,2,3),c(1,4,2,3)]
rownames(fst_matrix_ord) <- c("Chatham\nIsland", "Stewart\nIsland", "Fortrose", "Lincoln")
colnames(fst_matrix_ord) <- c("Chatham\nIsland", "Stewart\nIsland", "Fortrose", "Lincoln")

pheatmap(fst_matrix_ord, 
         color=viridis(50, option="A"), 
         na_col = "white",
         cluster_rows = FALSE,  # Disable row clustering
         cluster_cols = FALSE,
         border_color = NA,
         angle_col = 0,
         breaks = seq(0.016, 0.14, length.out = 51))

###########
## K = 2 ##
###########

# Read the PLINK2 fst matrix
fst_k2 <- fread("output/04_plink/fst/admix_k2/population_fst_results.fst.summary")


# Get unique population names
populations_k2 <- unique(c(fst_k2$`#POP1`, fst_k2$POP2))
# Initialize empty matrix
fst_matrix_k2 <- matrix(NA, nrow=length(populations_k2), ncol=length(populations_k2),
                     dimnames=list(populations_k2, populations_k2))

# Fill matrix using data from lower-triangle format
for (i in 1:nrow(fst_k2)) {
  pop1 <- fst_k2$`#POP1`[i]
  pop2 <- fst_k2$POP2[i]
  fst_value <- fst_k2$HUDSON_FST[i]
  
  fst_matrix_k2[pop1, pop2] <- fst_value
  fst_matrix_k2[pop2, pop1] <- fst_value  # Make symmetric
}

pheatmap(fst_matrix_k2, 
         color=viridis(50, option="A"), 
         na_col = "white",
         cluster_rows = FALSE,  # Disable row clustering
         cluster_cols = FALSE,
         border_color = NA,
         angle_col = 0,
         breaks = seq(0.016, 0.14, length.out = 51))

###########
## K = 3 ##
###########

# Read the PLINK2 fst matrix
fst_k3 <- fread("output/04_plink/fst/admix_k3/population_fst_results.fst.summary")


# Get unique population names
populations_k3 <- unique(c(fst_k3$`#POP1`, fst_k3$POP2))
# Initialize empty matrix
fst_matrix_k3 <- matrix(NA, nrow=length(populations_k3), ncol=length(populations_k3),
                        dimnames=list(populations_k3, populations_k3))

# Fill matrix using data from lower-triangle format
for (i in 1:nrow(fst_k3)) {
  pop1 <- fst_k3$`#POP1`[i]
  pop2 <- fst_k3$POP2[i]
  fst_value <- fst_k3$HUDSON_FST[i]
  
  fst_matrix_k3[pop1, pop2] <- fst_value
  fst_matrix_k3[pop2, pop1] <- fst_value  # Make symmetric
}

pheatmap(fst_matrix_k3, 
         color=viridis(50, option="A"), 
         na_col = "white",
         cluster_rows = FALSE,  # Disable row clustering
         cluster_cols = FALSE,
         border_color = NA,
         angle_col = 0,
         breaks = seq(0.016, 0.14, length.out = 51))

###########
## K = 4 ##
###########

# Read the PLINK2 fst matrix
fst_k4 <- fread("output/04_plink/fst/admix_k4/population_fst_results.fst.summary")


# Get unique population names
populations_k4 <- unique(c(fst_k4$`#POP1`, fst_k4$POP2))
# Initialize empty matrix
fst_matrix_k4 <- matrix(NA, nrow=length(populations_k4), ncol=length(populations_k4),
                        dimnames=list(populations_k4, populations_k4))

# Fill matrix using data from lower-triangle format
for (i in 1:nrow(fst_k4)) {
  pop1 <- fst_k4$`#POP1`[i]
  pop2 <- fst_k4$POP2[i]
  fst_value <- fst_k4$HUDSON_FST[i]
  
  fst_matrix_k4[pop1, pop2] <- fst_value
  fst_matrix_k4[pop2, pop1] <- fst_value  # Make symmetric
}

pheatmap(fst_matrix_k4, 
         color=viridis(50, option="A"), 
         na_col = "white",
         cluster_rows = FALSE,  # Disable row clustering
         cluster_cols = FALSE,
         border_color = NA,
         angle_col = 0,
         breaks = seq(0.016, 0.14, length.out = 51)) ## 50 colors → 51 breaks

