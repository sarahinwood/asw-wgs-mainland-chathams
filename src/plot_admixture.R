library(tidyverse)
library(ggplot2)
library(forcats)
library(data.table)
library(viridis)

## use clumpak aligned output for plotting so colours stay more consistent over different values of K

## make long table with all K res ###########

# your sample names
sample_names_table <- fread("output/04_plink/ld_pruned/ld_pruned_snps_plink_pca.fam", header = FALSE)
sample_names <- sample_names_table$V2

# list all Q files
q_files <- list.files("output/04_plink/ld_pruned/admixture/clumpak/output/1756952086/aligned.files", pattern = "admixture.*.Q.converted$", full.names = TRUE)

library(tidyverse)
library(stringr)

read_clumpak <- function(file, sample_names) {
  # extract K from filename (e.g. admixture.2.Q.converted → "2")
  k <- str_extract(basename(file), "(?<=\\.)\\d+(?=\\.Q)")
  k <- as.integer(k)
  
  # column names Q1..Qk
  qcols <- paste0("Q", seq_len(k))
  
  # read file as raw text, then split at ":" to keep only the right-hand side
  df <- read_delim(file, 
                   delim = ":", 
                   col_names = c("meta", "qvals"),
                   trim_ws = TRUE,
                   show_col_types = FALSE)
  
  # split qvals string into numeric columns
  df_q <- df %>%
    select(qvals) %>%
    separate(qvals, into = qcols, sep = "\\s+", convert = TRUE, extra = "drop")
  
  # add metadata
  df_q$sample <- sample_names
  df_q$k <- k
  
  # pivot to long format
  df_long <- df_q %>%
    pivot_longer(cols = all_of(qcols), names_to = "Q", values_to = "value")
  
  # assign Location based on sample name
  df_long$Location <- case_when(
    grepl("Chatham", df_long$sample) ~ "Chatham Island",
    grepl("Stewart", df_long$sample) ~ "Stewart Island",
    grepl("Lincoln", df_long$sample) ~ "Lincoln",
    grepl("Fortrose", df_long$sample) ~ "Fortrose",
    TRUE ~ NA_character_
  )
  
  return(df_long)
}

# apply to all files and combine
long_admix_all <- map_dfr(q_files, ~read_clumpak(.x, sample_names))

# Fix sample order across all Ks (by Location then sample name)
sample_order <- long_admix_all %>%
  distinct(sample, Location) %>%
  arrange(Location, sample) %>%
  pull(sample)
long_admix_all$sample <- factor(long_admix_all$sample, levels = sample_order)
long_admix_all$k <- paste("K =", long_admix_all$k)
long_admix_all$k <- factor(long_admix_all$k, levels=c("K = 2",
                                                      "K = 3",
                                                      "K = 4",
                                                      "K = 5",
                                                      "K = 6",
                                                      "K = 7",
                                                      "K = 8",
                                                      "K = 9",
                                                      "K = 10"))

long_admix_all$Location <- factor(long_admix_all$Location, levels=c("Chatham Island",
                                                                    "Stewart Island",
                                                                    "Fortrose",
                                                                    "Lincoln"))

#fwrite(long_admix_all, "output/04_plink/ld_pruned/admixture/clumpak/output/admix_clumpak_table_for_plot.csv")

## plot all values of K at once
# coloured manually to match manual pca
manual_colours = c("Q1"="#0d0887",
                   "Q2"="#f0f921",
                   "Q3"="#bd3786",
                   "Q4"="#fb9f3a",
                   "Q5"="#7201a8",
                   "Q6"="#fdca26",
                   "Q7"="#ed7953",
                   "Q8"="#46039f",
                   "Q9"="#d8576b",
                   "Q10"="#9c179e")


ggplot(long_admix_all, aes(x = sample, y = value, fill = factor(Q))) + 
  geom_col(width = 1) +
  facet_grid(k ~ Location, scales = "free_x", space = "free_x", switch = "x") +
  xlab("Population") +
  ylab("Ancestry") +
  scale_y_continuous(limits = c(0, 1.0001),
                     breaks = c(0, 0.5, 1.0),
                     labels = c("0", "0.5", "1.0")) +
  theme_minimal() +
  scale_fill_manual(values=manual_colours) +
  theme(
    axis.ticks.y = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing.x = unit(1, "mm"),
    axis.text.x = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(size = 10, face = "bold"),
    strip.text.y = element_text(angle = 0),   # <--- makes k horizontal
    strip.switch.pad.grid = unit(0.001, "mm"),
    plot.margin = margin(5, 5, 5, 5),
    legend.position = "none"
  )

#### sorting samples by Q proportions while keeping same order throughout plots #### 

# K=2: Q1 and Q2
df_K2 <- long_admix_all %>%
  filter(k == "K = 2") %>%
  select(sample, Location, Q, value) %>%
  pivot_wider(names_from = Q, values_from = value) %>%
  rename(Q1_K2 = Q1, Q2_K2 = Q2)

# K=3: Q2
df_K3 <- long_admix_all %>%
  filter(k == "K = 3", Q == "Q2") %>%
  select(sample, Location, value) %>%
  rename(Q2_K3 = value)

# K=4: Q4
df_K4 <- long_admix_all %>%
  filter(k == "K = 4", Q == "Q1") %>%
  select(sample, Location, value) %>%
  rename(Q4_K4 = value)

# K=7: Q1
df_K7 <- long_admix_all %>%
  filter(k == "K = 7", Q == "Q1") %>%
  select(sample, Location, value) %>%
  rename(Q1_K7 = value)

# K=9: Q9
df_K9 <- long_admix_all %>%
  filter(k == "K = 9", Q == "Q9") %>%
  select(sample, Location, value) %>%
  rename(Q9_K9 = value)

# K=10: Q9
df_K10 <- long_admix_all %>%
  filter(k == "K = 10", Q == "Q9") %>%
  select(sample, Location, value) %>%
  rename(Q1_K10 = value)

# Join all together
sort_df <- df_K2 %>%
  left_join(df_K3, by = c("sample", "Location")) %>%
  left_join(df_K4, by = c("sample", "Location")) %>%
  left_join(df_K7, by = c("sample", "Location")) %>%
  left_join(df_K9, by = c("sample", "Location")) %>%
  left_join(df_K10, by = c("sample", "Location")) %>%
  mutate(
    Q2_K3   = ifelse(is.na(Q2_K3), 0, Q2_K3),
    Q4_K4   = ifelse(is.na(Q4_K4), 0, Q4_K4),
    Q1_K7   = ifelse(is.na(Q1_K7), 0, Q1_K7),
    Q9_K9   = ifelse(is.na(Q9_K9), 0, Q9_K9),
    Q1_K10  = ifelse(is.na(Q1_K10), 0, Q1_K10)
  )

# Assign large constants to ensure deterministic priority
sort_df <- sort_df %>%
  mutate(sort_key = Q1_K2*1e15 +   # highest priority
           Q2_K2*1e12 +
           Q2_K3*1e9 +
           Q4_K4*1e6 +
           Q1_K7*1e3 +
           Q9_K9*1 +
           Q1_K10*1e-3)   # lowest priority tie-breaker
sample_order <- sort_df %>%
  group_by(Location) %>%
  arrange(desc(sort_key), .by_group = TRUE) %>%
  pull(sample)

long_admix_sorted <- long_admix_all %>%
  mutate(sample = factor(sample, levels = sample_order))

long_admix_sorted$Q <- factor(long_admix_sorted$Q, levels=c("Q1",
                                                            "Q2",
                                                            "Q3",
                                                            "Q4",
                                                            "Q5",
                                                            "Q6",
                                                            "Q7",
                                                            "Q8",
                                                            "Q9",
                                                            "Q10"))

ggplot(long_admix_sorted, aes(x = sample, y = value, fill = Q)) +
  geom_col(width = 1) +
  facet_grid(k ~ Location, scales = "free_x", space = "free_x", switch = "x") +
  xlab("Population") +
  ylab("") +
  scale_y_continuous(limits = c(0, 1.0001),
                     breaks = c(0, 0.5, 1.0),
                     labels = c("0", "0.5", "1.0")) +
  theme_minimal() +
  scale_fill_manual(values = manual_colours) +
  theme(
    axis.ticks.y = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing.x = unit(1, "mm"),
    axis.text.x = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(size = 10, face = "bold"),
    strip.text.y = element_text(angle = 0),   # <--- makes k horizontal
    strip.switch.pad.grid = unit(0.001, "mm"),
    plot.margin = margin(5, 5, 5, 5),
    legend.position = "none"
  )


subset_k <- c("K = 2", "K = 3", "K = 4")
long_admix_sorted_k234 <- subset(long_admix_sorted, k %in% subset_k)
ggplot(long_admix_sorted_k234, aes(x = sample, y = value, fill = Q)) +
  geom_col(width = 1) +
  facet_grid(k ~ Location, scales = "free_x", space = "free_x", switch = "x") +
  xlab("Population") +
  ylab("") +
  scale_y_continuous(limits = c(0, 1.0001),
                     breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                     labels = c("0", "0.25", "0.5", "0.75", "1.0")) +
  theme_minimal() +
  scale_fill_manual(values = manual_colours) +
  theme(
    axis.ticks.y = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing.x = unit(1, "mm"),
    axis.text.x = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(size = 10, face = "bold"),
    strip.text.y = element_text(angle = 0),   # <--- makes k horizontal
    strip.switch.pad.grid = unit(0.001, "mm"),
    plot.margin = margin(5, 5, 5, 5),
    legend.position = "none"
  )










## sort bars on proportion so highest at bottom
df_stacked <- long_admix_sorted %>%
  group_by(k, Location, sample) %>%
  arrange(desc(value), .by_group = TRUE) %>%  # largest value first = bottom
  mutate(
    ymin = c(0, head(cumsum(value), -1)),
    ymax = cumsum(value),
    xi = as.integer(sample)  # numeric x for geom_rect
  ) %>%
  ungroup()

ggplot(df_stacked, aes(ymin = ymin, ymax = ymax, xmin = xi - 0.5, xmax = xi + 0.5, fill = Q)) +
  geom_rect() +
  scale_x_continuous(
    breaks = unique(df_stacked$xi),
    labels = levels(df_stacked$sample),
    expand = c(0,0)
  ) +
  facet_grid(k ~ Location, scales = "free_x", space = "free_x", switch = "x") +
  xlab("Population") +
  ylab("Ancestry") +
  scale_y_continuous(limits = c(0, 1.0001), breaks = c(0, 0.5, 1.0), labels = c("0","0.5","1.0")) +
  theme_minimal() +
  scale_fill_manual(values = manual_colours) +
  theme(
    axis.text.x = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(size = 10, face = "bold"),
    strip.text.y = element_text(angle = 0),
    panel.grid = element_blank(),
    legend.position = "none"
  )
