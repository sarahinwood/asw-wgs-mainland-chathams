library(data.table)
library(tidyverse)

##############
## tajima d ##
##############

## read in Tajimas D output
tajd_chat <- fread('output/06_pop_stats/tajimas_d_Chatham_100kb.Tajima.D', na.strings="NaN")
tajd_chat$Location <- paste("Chatham Island")

tajd_stew <- fread('output/06_pop_stats/tajimas_d_Stewart_100kb.Tajima.D', na.strings="NaN")
tajd_stew$Location <- paste("Stewart Island")

tajd_fort <- fread('output/06_pop_stats/tajimas_d_Fortrose_100kb.Tajima.D', na.strings="NaN")
tajd_fort$Location <- paste("Fortrose")

tajd_linc <- fread('output/06_pop_stats/tajimas_d_Lincoln_100kb.Tajima.D', na.strings="NaN")
tajd_linc$Location <- paste("Lincoln")

## join
tajd_full <- full_join(tajd_chat, full_join(tajd_stew, full_join(tajd_fort, tajd_linc)))
tajd_full[] <- lapply(tajd_full, function(x) if (is.numeric(x)) replace(x, is.nan(x), NA) else x)

## reorder Location as factor
tajd_full$Location <- factor(tajd_full$Location, 
                           levels=c("Chatham Island", "Stewart Island", "Fortrose", "Lincoln"))

ggplot(tajd_full, aes(x=Location, y=TajimaD, fill=Location))+
  geom_violin()+
  scale_fill_viridis(discrete=T, alpha=0.9)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5)+
  theme_minimal()

tajd_full %>%
  group_by(Location) %>%
  summarise(mean=mean(TajimaD, na.rm = TRUE),
            median=median(TajimaD, na.rm = TRUE),
            min=min(TajimaD, na.rm = TRUE),
            max=max(TajimaD, na.rm = TRUE))
# positive Tajima's D signifies low levels of both low and high frequency polymorphisms/excess of common variation in a Location, indicating a decrease in population size and/or balancing selection
# Negative values of Tajima's D, on the other hand, indicate an excess of rare variation, consistent with population growth, or positive selection
