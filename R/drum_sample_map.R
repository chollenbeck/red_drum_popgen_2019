library(tidyverse)
library(ggmap)
library(adegenet)

coords <- read.csv("data/raw/sample_coords.csv")

# Read in the population data
gen <- read.genepop("data/raw/genepop.18.gen", ncode=3)

# Assign the individuals names to the genind object
num.loci <- length(locNames(gen))
ind.tab <- read.table("data/raw/genepop.18.gen", skip=num.loci+2, fill=TRUE)
ind.tab <- subset(ind.tab, V1 !='POP')
inds <- ind.tab$V1
inds <- gsub(',', '', inds)
indNames(gen) <- inds

# Read in the library data
lib <- read.table("data/raw/drum_sample_data.txt", header=TRUE, sep="\t")
#samp <- read.table("data/external/RD_Master_Sample_List.txt", header=TRUE, fill=TRUE, sep="\t")

samples <- lib %>%
                filter(Sample %in% inds) %>%
                group_by(Locality) %>%
                summarise(n.inds = n())

colnames(samples)[1] <- "name"

# Filter out the replicate samples from the table
samples <-samples %>% filter(name != "LLM" & name != "MAT")
samples$n.inds[6] <- "22/34"
samples$n.inds[7] <- "40/48"
samples$name[6] <- "LLM"
samples$name[7] <- "MAT"



labels <- c("LLM", "MAT", "SAB", "MIS", "APA", "CEK", "CHA", "IND", "HAR", "WAS", "SCA")

#samples$name <- factor(samples$name, levels = labels)

sample_coords <- left_join(samples, coords, by = "name")


states <- map_data("state")
mexico <- map_data("world") %>% filter(region == "Mexico")

usa_mex <- rbind(states, mexico)

state_lab <- tibble::tribble(
  ~state, ~x, ~y,
  "TX", -96, 31,
  "LA", -92.5, 32,
  "MS", -89.5, 33,
  "AL", -86.5, 33,
  "FL", -81.5, 28,
  "GA", -83, 32.5,
  "SC", -81, 33.5)

ggsave(filename = "out/fig/figure_1.png", plot = map_plot, width = 24, height = 16, units = "cm") 


