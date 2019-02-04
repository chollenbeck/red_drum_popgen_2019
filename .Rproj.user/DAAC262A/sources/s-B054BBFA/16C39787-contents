library(dplyr)
library(ggplot2)
library(adegenet)


all_dat <- read.table("data/raw/full_outlier_summary.txt", header = TRUE)


# Plot FST by heterozygosity, with points colored by number of methods of detection

het_fst_plot <- ggplot(data=all_dat, aes(x=het, y=fst)) +
  geom_point(aes(colour=as.factor(num.meth)), size=4, alpha = 0.8) +
  geom_hline(yintercept=0, colour="red") +
  scale_colour_manual(name="Number of Methods",
                      values=c("#848484", "#00441b", "#4575b4", "#b2182b")) +
  labs(x = "Heterozygosity", y = expression(F[ST])) +
  rd_theme

ggsave(filename = "out/fig/figure_S2.png", plot = het_fst_plot, width = 36, height = 24, units = "cm")

