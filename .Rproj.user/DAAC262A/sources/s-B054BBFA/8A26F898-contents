#### Produce a multi-panel figure summarizing the population structure


library(adegenet)
library(ggplot2)
library(dplyr)
library(gridExtra)


### Shared variables

pca_col <- c('#9e0142','#d53e4f','#f46d43','#fdae61','#a1d99b','#41ab5d','#00441b','#4eb3d3','#2b8cbe','#0868ac','#084081')


### PCA - neutral only ---------------------------

neut_gen <- read.genepop("data/raw/final_neutral.gen", ncode = 3)


# Plot a PCA using only neutral loci
x_neut <- scaleGen(neut_gen, NA.method="mean")
neut_pca <- dudi.pca(x_neut,center=FALSE,scale=FALSE,scannf=FALSE,nf=4)

eig_percent_neut <- round((neut_pca$eig/(sum(neut_pca$eig)))*100,2)

pca_dat_neut <- as.data.frame(neut_pca$li)
pca_dat_neut <- cbind(pca_dat_neut, pop = pop(neut_gen))
pca_dat_neut$pop <- gsub("_\\d+", "", pca_dat_neut$pop, perl=TRUE)
pca_dat_neut$pop <- gsub("LM2", "LLM", pca_dat_neut$pop, perl=TRUE)
pca_dat_neut$pop <- gsub("MB2", "MAT", pca_dat_neut$pop, perl=TRUE)
pca_dat_neut$pop <- factor(pca_dat_neut$pop, levels = c("LLM", "MAT", "SAB", "MIS", "APA", "CEK", "CHA", "IND", "HAR", "WAS", "SCA"))

pca_dat_neut <- mutate(pca_dat_neut, region = case_when(pop %in% c("LLM", "MAT", "SAB", "MIS") ~ "NWG",
                       pop %in% c("APA", "CEK", "CHA") ~ "NEG",
                       pop %in% c("IND", "HAR", "WAS", "SCA") ~ "ATL")) %>%
  mutate(region = factor(region, levels = c("NWG", "NEG", "ATL")))

# Get just the PCA legend
pca_legend <- get_legend(
    ggplot(data=pca_dat_neut) +
    geom_point(aes(x=Axis1, y=Axis2, colour=pop, shape = region), size=12, alpha=0.8) +
    scale_colour_manual(values = pca_col, name = "Locality") +
      scale_shape_manual(breaks = c("NWG", "NEG", "ATL"), values = c(16, 17, 18), name = "Region") +
    rd_theme +
    theme(axis.title.y = element_text(vjust = 0.5, angle = 90),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 20, face = "bold")) +
    labs(x = paste("", "Axis 1: ", eig_percent_neut[1], "%"), y = paste("", "Axis 2: ", eig_percent_neut[2], "%"))
)



pca_neut <- ggplot(data=pca_dat_neut) +
  geom_point(aes(x=Axis1, y=Axis2, colour=pop, shape = region), size=7, alpha=0.8) +
  scale_colour_manual(values = pca_col) +
  scale_shape_manual(breaks = c("NWG", "NEG", "ATL"), values = c(16, 17, 18)) +
  rd_theme +
  theme(legend.position = "none", axis.title.y = element_text(vjust = 0.5, angle = 90)) +
  labs(x = paste("", "Axis 1: ", eig_percent_neut[1], "%"), y = paste("", "Axis 2: ", eig_percent_neut[2], "%")) +
  annotate(geom = "text", label = "a)", x = -13, y = 16, size = 15)



### PCA - outlier only ---------------------------

out_gen <- read.genepop("data/raw/final_outliers.gen", ncode = 3)


# Plot a PCA using only outlier loci
x_out <- scaleGen(out_gen, NA.method="mean")
out_pca <- dudi.pca(x_out,center=FALSE,scale=FALSE,scannf=FALSE,nf=4)

eig_percent_out <- round((out_pca$eig/(sum(out_pca$eig)))*100,2)

pca_dat_out <- as.data.frame(out_pca$li)
pca_dat_out <- cbind(pca_dat_out, pop = pop(out_gen))
pca_dat_out$pop <- gsub("_\\d+", "", pca_dat_out$pop, perl=TRUE)
pca_dat_out$pop <- gsub("LM2", "LLM", pca_dat_out$pop, perl=TRUE)
pca_dat_out$pop <- gsub("MB2", "MAT", pca_dat_out$pop, perl=TRUE)
pca_dat_out$pop <- factor(pca_dat_out$pop, levels = c("LLM", "MAT", "SAB", "MIS", "APA", "CEK", "CHA", "IND", "HAR", "WAS", "SCA"))

pca_dat_out <- mutate(pca_dat_out, region = case_when(pop %in% c("LLM", "MAT", "SAB", "MIS") ~ "NWG",
                                                        pop %in% c("APA", "CEK", "CHA") ~ "NEG",
                                                        pop %in% c("IND", "HAR", "WAS", "SCA") ~ "ATL")) %>%
  mutate(region = factor(region, levels = c("NWG", "NEG", "ATL")))

pca_out <- ggplot(data=pca_dat_out) +
  geom_point(aes(x=Axis1, y=Axis2*-1, colour=pop, shape = region), size=7, alpha=0.8) +
  scale_colour_manual(values = pca_col) +
  scale_shape_manual(breaks = c("NWG", "NEG", "ATL"), values = c(16, 17, 18)) +
  rd_theme +
  theme(legend.position = "none", axis.title.y = element_text(vjust = 0.5, angle = 90)) +
  labs(x = paste("", "Axis 1: ", eig_percent_out[1], "%"), y = paste("", "Axis 2: ", eig_percent_out[2], "%")) +
  annotate(geom = "text", label = "d)", x = -8, y = 7, size = 15)

### Pairwise Fst - neutral loci ---------------------------

# Read in the pairwise Fst data

neut_fst_tab <- read.table("data/raw/neut_fst_tbl.txt", header = TRUE, sep = "\t")

# Order the factors geographically for plotting the data

neut_fst_tab$site1 <- factor(neut_fst_tab$site1, c("LLM", "MAT", "SAB", "MIS", "APA", "CEK", "CHA", "IND", "HAR", "WAS", "SCA"))
neut_fst_tab$site2 <- factor(neut_fst_tab$site2, c("LLM", "MAT", "SAB", "MIS", "APA", "CEK", "CHA", "IND", "HAR", "WAS", "SCA"))


neut_fst_tab <- mutate(neut_fst_tab, "region.1.wgulf" = (site1 == "LLM" | site1 == "MAT" | site1 == "SAB" | site1 == "MIS"))
neut_fst_tab <- mutate(neut_fst_tab, "region.2.wgulf" = (site2 == "LLM" | site2 == "MAT" | site2 == "SAB" | site2 == "MIS"))
neut_fst_tab <- mutate(neut_fst_tab, "region.1.egulf" = (site1 == "APA" | site1 == "CEK" | site1 == "CHA"))
neut_fst_tab <- mutate(neut_fst_tab, "region.2.egulf" = (site2 == "APA" | site2 == "CEK" | site2 == "CHA"))
neut_fst_tab <- mutate(neut_fst_tab, "region.1.atl" = (site1 == "IND" | site1 == "HAR" | site1 == "WAS" | site1 == "SCA"))
neut_fst_tab <- mutate(neut_fst_tab, "region.2.atl" = (site2 == "IND" | site2 == "HAR" | site2 == "WAS" | site2 == "SCA"))

neut_fst_tab <- mutate(neut_fst_tab, "wgulf.egulf" = ((region.1.wgulf == TRUE & region.2.egulf == TRUE)) | (region.2.wgulf == TRUE & region.1.egulf == TRUE))
neut_fst_tab <- mutate(neut_fst_tab, "wgulf.atl" = (region.1.wgulf == TRUE & region.2.atl == TRUE) | (region.2.wgulf == TRUE & region.1.atl == TRUE))
neut_fst_tab <- mutate(neut_fst_tab, "egulf.atl" = (region.1.egulf == TRUE & region.2.atl == TRUE) | (region.2.egulf == TRUE & region.1.atl == TRUE))

neut_fst_tab$comparison[neut_fst_tab$wgulf.egulf == TRUE] <- "WGulf.EGulf"
neut_fst_tab$comparison[neut_fst_tab$wgulf.atl == TRUE] <- "WGulf.Atl"
neut_fst_tab$comparison[neut_fst_tab$egulf.atl == TRUE] <- "EGulf.Atl"
neut_fst_tab$comparison[neut_fst_tab$region.1.wgulf == TRUE & neut_fst_tab$region.2.wgulf == TRUE] <- "Within Region"
neut_fst_tab$comparison[neut_fst_tab$region.1.egulf == TRUE & neut_fst_tab$region.2.egulf == TRUE] <- "Within Region"
neut_fst_tab$comparison[neut_fst_tab$region.1.atl == TRUE & neut_fst_tab$region.2.atl == TRUE] <- "Within Region"



# Read in the Mantel test data
mantel <- read.table("data/raw/mantel.txt", header = TRUE, sep = "\t")
mant_stat_neut <- mantel$Statistic[1]
mant_pval_neut <- mantel$P.value[1]


# Plot an FST heatmap

heatmap_neut <- ggplot(neut_fst_tab, aes(site1, site2)) +
  geom_tile(aes(fill=fst)) +
  scale_fill_continuous(low="blue", high="red", name=expression(F[ST])) +
  rd_theme +
  theme(axis.text=element_text(face="bold"),
        legend.text = element_text(color="black", size=16, angle = 0),
        legend.title = element_text(color="black", size=20, angle = 0),
        axis.text.x = element_text(angle = 90),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
    annotate(geom = "text", label = "b)", x = 1.25, y = 10.75, size = 15)


# Get just the legend for the IBD figure
ibd_legend <- get_legend(
    ggplot(neut_fst_tab, aes(x=dist, y=fst)) +
        geom_point(aes(colour=comparison, shape = comparison), size=14, alpha=0.8) +
        geom_smooth(method=lm, se=FALSE, colour="black", size=1.5) +
        #geom_smooth(data=filter(fst.tab, comparison=="EGulf.Atl"), aes(x=dist, y=fst), method=lm, se=FALSE, colour="#00441b") +
        #geom_smooth(data=filter(fst.tab, comparison=="WGulf.EGulf"), aes(x=dist, y=fst), method=lm, se=FALSE, colour="#9e0142") +
        #geom_smooth(data=filter(fst.tab, comparison=="WGulf.Atl"), aes(x=dist, y=fst), method=lm, se=FALSE, colour="#084081") +
        scale_colour_manual(name="Comparison",
                            breaks=c("Within Region", "EGulf.Atl", "WGulf.EGulf", "WGulf.Atl"),
                            labels=c("Within Region", "NEG-ATL", "NWG-NEG", "NWG-ATL"),
                            values=c("#00441b", "#084081", "#9e0142", "darkgrey")) +
        scale_shape_manual(name="Comparison",
                         breaks=c("Within Region", "EGulf.Atl", "WGulf.EGulf", "WGulf.Atl"),
                         labels=c("Within Region", "NEG-ATL", "NWG-NEG", "NWG-ATL"),
                         values=c(15, 16, 17, 18)) +
        rd_theme +
        labs(x = "Geographic distance (km)", y = expression(F[ST])) +
        #geom_text(aes(x=3000, y=0.007, label=paste("Mantel r: ", sprintf("%.3f", mant$statistic), " (P-value = ", sprintf("%.3f", mant$signif), ")", sep=''))) +
        theme(axis.title.y = element_text(vjust = 0.5, angle = 90),
              legend.text = element_text(size = 18),
              legend.title = element_text(size = 20, face = "bold"))
)

# Plot Fst against geographic distance, coloring the points by regional comparison
ibd_neut <- ggplot(neut_fst_tab, aes(x=dist, y=fst)) +
  geom_point(aes(colour=comparison, shape=comparison), size=10, alpha=0.8) +
  geom_smooth(method=lm, se=FALSE, colour="black", size=1.5) +
  #geom_smooth(data=filter(fst.tab, comparison=="EGulf.Atl"), aes(x=dist, y=fst), method=lm, se=FALSE, colour="#00441b") +
  #geom_smooth(data=filter(fst.tab, comparison=="WGulf.EGulf"), aes(x=dist, y=fst), method=lm, se=FALSE, colour="#9e0142") +
  #geom_smooth(data=filter(fst.tab, comparison=="WGulf.Atl"), aes(x=dist, y=fst), method=lm, se=FALSE, colour="#084081") +
  scale_colour_manual(name="Comparison",
                      breaks=c("Within Region", "EGulf.Atl", "WGulf.EGulf", "WGulf.Atl"),
                      labels=c("Within Region", "NEG-ATL", "NWG-NEG", "NWG-ATL"),
                      values=c("#00441b", "#084081", "#9e0142", "darkgrey")) +
  scale_shape_manual(name="Comparison",
                      breaks=c("Within Region", "EGulf.Atl", "WGulf.EGulf", "WGulf.Atl"),
                      labels=c("Within Region", "NEG-ATL", "NWG-NEG", "NWG-ATL"),
                      values=c(15, 16, 17, 18)) +
  rd_theme +
  labs(x = "Geographic distance (km)", y = expression(F[ST])) +
  #geom_text(aes(x=3000, y=0.007, label=paste("Mantel r: ", sprintf("%.3f", mant$statistic), " (P-value = ", sprintf("%.3f", mant$signif), ")", sep=''))) +
  theme(legend.position="none", axis.title.y = element_text(vjust = 0.5, angle = 90)) +
  annotate(geom = "text", label = "c)", x = 250, y = 0.0068, size = 15)



### Pairwise Fst - outlier loci ---------------------------

# Read in the pairwise Fst data

out_fst_tab <- read.table("data/raw/out_fst_tbl.txt", header = TRUE, sep = "\t")


# Order the factors geographically for plotting the data

out_fst_tab$site1 <- factor(out_fst_tab$site1, c("LLM", "MAT", "SAB", "MIS", "APA", "CEK", "CHA", "IND", "HAR", "WAS", "SCA"))
out_fst_tab$site2 <- factor(out_fst_tab$site2, c("LLM", "MAT", "SAB", "MIS", "APA", "CEK", "CHA", "IND", "HAR", "WAS", "SCA"))


out_fst_tab <- mutate(out_fst_tab, "region.1.wgulf" = (site1 == "LLM" | site1 == "MAT" | site1 == "SAB" | site1 == "MIS"))
out_fst_tab <- mutate(out_fst_tab, "region.2.wgulf" = (site2 == "LLM" | site2 == "MAT" | site2 == "SAB" | site2 == "MIS"))
out_fst_tab <- mutate(out_fst_tab, "region.1.egulf" = (site1 == "APA" | site1 == "CEK" | site1 == "CHA"))
out_fst_tab <- mutate(out_fst_tab, "region.2.egulf" = (site2 == "APA" | site2 == "CEK" | site2 == "CHA"))
out_fst_tab <- mutate(out_fst_tab, "region.1.atl" = (site1 == "IND" | site1 == "HAR" | site1 == "WAS" | site1 == "SCA"))
out_fst_tab <- mutate(out_fst_tab, "region.2.atl" = (site2 == "IND" | site2 == "HAR" | site2 == "WAS" | site2 == "SCA"))

out_fst_tab <- mutate(out_fst_tab, "wgulf.egulf" = ((region.1.wgulf == TRUE & region.2.egulf == TRUE)) | (region.2.wgulf == TRUE & region.1.egulf == TRUE))
out_fst_tab <- mutate(out_fst_tab, "wgulf.atl" = (region.1.wgulf == TRUE & region.2.atl == TRUE) | (region.2.wgulf == TRUE & region.1.atl == TRUE))
out_fst_tab <- mutate(out_fst_tab, "egulf.atl" = (region.1.egulf == TRUE & region.2.atl == TRUE) | (region.2.egulf == TRUE & region.1.atl == TRUE))

out_fst_tab$comparison[out_fst_tab$wgulf.egulf == TRUE] <- "WGulf.EGulf"
out_fst_tab$comparison[out_fst_tab$wgulf.atl == TRUE] <- "WGulf.Atl"
out_fst_tab$comparison[out_fst_tab$egulf.atl == TRUE] <- "EGulf.Atl"
out_fst_tab$comparison[out_fst_tab$region.1.wgulf == TRUE & out_fst_tab$region.2.wgulf == TRUE] <- "Within Region"
out_fst_tab$comparison[out_fst_tab$region.1.egulf == TRUE & out_fst_tab$region.2.egulf == TRUE] <- "Within Region"
out_fst_tab$comparison[out_fst_tab$region.1.atl == TRUE & out_fst_tab$region.2.atl == TRUE] <- "Within Region"

# Read in the Mantel test data
mantel <- read.table("data/raw/mantel.txt", header = TRUE, sep = "\t")
mant_stat_out <- mantel$Statistic[2]
mant_pval_out <- mantel$P.value[2]


# Plot an FST heatmap

heatmap_out <- ggplot(out_fst_tab, aes(site1, site2)) +
    geom_tile(aes(fill=fst)) +
    scale_fill_continuous(low="blue", high="red", name=expression(F[ST])) +
    rd_theme +
    theme(axis.text=element_text(face="bold"),
          legend.text = element_text(color="black", size=16, angle = 0),
          legend.title = element_text(color="black", size=20, angle = 0),
          axis.text.x = element_text(angle = 90),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    annotate(geom = "text", label = "e)", x = 1.25, y = 10.75, size = 15)




# Plot Fst against geographic distance, coloring the points by regional comparison
ibd_out <- ggplot(out_fst_tab, aes(x=dist, y=fst)) +
  geom_point(aes(colour = comparison, shape = comparison), size=10, alpha=0.8) +
  geom_smooth(method=lm, se=FALSE, colour="black", size=1.5) +
  #geom_smooth(data=filter(fst.tab, comparison=="EGulf.Atl"), aes(x=dist, y=fst), method=lm, se=FALSE, colour="#00441b") +
  #geom_smooth(data=filter(fst.tab, comparison=="WGulf.EGulf"), aes(x=dist, y=fst), method=lm, se=FALSE, colour="#9e0142") +
  #geom_smooth(data=filter(fst.tab, comparison=="WGulf.Atl"), aes(x=dist, y=fst), method=lm, se=FALSE, colour="#084081") +
  scale_colour_manual(name="Comparison",
                      breaks=c("Within Region", "EGulf.Atl", "WGulf.EGulf", "WGulf.Atl"),
                      labels=c("Within Region", "NEG-ATL", "NWG-NEG", "NWG-ATL"),
                      values=c("#00441b", "#084081", "#9e0142", "darkgrey")) +
  scale_shape_manual(name="Comparison",
                     breaks=c("Within Region", "EGulf.Atl", "WGulf.EGulf", "WGulf.Atl"),
                     labels=c("Within Region", "NEG-ATL", "NWG-NEG", "NWG-ATL"),
                     values=c(15, 16, 17, 18)) +
  rd_theme +
  labs(x = "Geographic distance (km)", y = expression(F[ST])) +
  #geom_text(aes(x=3000, y=0.14, label=paste("Mantel r: ", sprintf("%.3f", mant$statistic), " (P-value = ", sprintf("%.3f", mant$signif), ")", sep=''))) +
  theme(legend.position="none", axis.title.y = element_text(vjust = 0.5, angle = 90)) +
  annotate(geom = "text", label = "f)", x = 250, y = 0.23, size = 15)


grid.arrange(pca_legend, pca_neut, pca_out, heatmap_neut, heatmap_out, ibd_neut, ibd_out, ibd_legend,
             layout_matrix = rbind(c(1,2,4,6,8), c(1,3,5,7,8)),
             ncol=5, widths = c(5, 10, 10, 10, 5), heights = c(12.5, 12.5))


### Make the combined plot

pop_str_plot <- grid.arrange(pca_legend, pca_neut, pca_out, heatmap_neut, heatmap_out, ibd_neut, ibd_out, ibd_legend,
             layout_matrix = rbind(c(1,2,4,6,8), c(1,3,5,7,8)),
             ncol=5, widths = c(2, 12, 15, 12, 2), heights = c(12.5, 12.5))



ggsave(filename = "out/fig/figure_2.png", plot = pop_str_plot, width = 64, height = 28, units = "cm") 

