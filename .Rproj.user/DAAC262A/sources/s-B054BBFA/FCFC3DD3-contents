library(tidyverse)
library(adegenet)
library(hierfstat)

gen <- read.genepop("data/raw/final_outliers.gen", ncode = 3)

# Fix pop names
popNames(gen) <- popNames(gen) %>%
  str_replace("_\\d+", "") %>%
  str_replace("LM2", "LLM") %>%
  str_replace("MB2", "MAT")


# Code by region
reg <- gen

# Recode the pops by region
popNames(reg) <- popNames(reg) %>%
  as.data.frame() %>%
  mutate(pop = case_when(`.` %in% c("LLM", "MAT", "SAB", "MIS") ~ "NWG",
                         `.` %in% c("APA", "CEK", "CHA") ~ "NEG",
                         `.` %in% c("IND", "HAR", "WAS", "SCA") ~ "ATL")) %>%
  pull(pop)

# Get allele freqs
mca_plot <- reg %>%
  genind2genpop() %>%
  makefreq() %>%
  as.data.frame() %>%
  rownames_to_column(var = "pop") %>%
  gather(locus.allele, freq, -pop) %>%
  separate(locus.allele, into = c("locus", "allele"), sep = "\\.") %>%
  group_by(pop, locus) %>%
  slice(which.max(freq)) %>%
  ggplot(aes(x = freq, fill = pop)) +
    geom_density(alpha = 0.5, adjust = 1.5) +
    theme_minimal() +
    labs(x = "Frequency of most common allele (MCA) by region",
         fill = "Region") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank())

ggsave(filename = "out/fig/figure_3.png", plot = mca_plot, width = 24, height = 18, units = "cm") 


