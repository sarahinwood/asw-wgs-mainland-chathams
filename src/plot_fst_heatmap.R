library(data.table)
library(tidyverse)
library(reshape2)
library(pheatmap)
library(viridis)
library(reshape2)

# Read the PLINK2 fst matrix
fst <- fread("output/04_plink/fst/population_fst_results.fst.summary")

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
         angle_col = 0)