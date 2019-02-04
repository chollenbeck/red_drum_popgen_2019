library(adegenet)
library(ggplot2)
library(dplyr)

### ggplot PCA - including temporal sample ---------------------------

gen <- read.genepop("data/raw/genepop.18.gen", ncode=3)

# Assign the individuals names to the genind object
num.loci <- length(locNames(gen))
ind.tab <- read.table("data/raw/genepop.18.gen", skip=num.loci+2, fill=TRUE)
ind.tab <- subset(ind.tab, V1 !='POP')
inds <- ind.tab$V1
inds <- gsub(',', '', inds)
indNames(gen) <- inds


# Generate a PCA using all of the samples with adegenet

x <- scaleGen(gen, NA.method="mean")
pca <- dudi.pca(x,center=FALSE,scale=FALSE,scannf=FALSE,nf=4)

eig_percent <- round((pca$eig/(sum(pca$eig)))*100,2)


# Convert the results into a dataframe compatible with ggplot2

pca.dat <- as.data.frame(pca$li)
pca.dat <- cbind(pca.dat, pop = pop(gen))
pca.dat$pop <- gsub("_\\d+", "", pca.dat$pop, perl=TRUE)
pca.dat$pop <- factor(pca.dat$pop, levels = c("LLM", "LM2", "MAT", "MB2", "SAB", "MIS", "APA", "CEK", "CHA", "IND", "HAR", "WAS", "SCA"))

col <- c('#9e0142', '#9e0142','#d53e4f', '#d53e4f', '#f46d43','#fdae61','#a1d99b','#41ab5d','#00441b','#4eb3d3','#2b8cbe','#0868ac','#084081')
labels <- c("LLM (2008)", "LLM (2014/2015)", "MAT (2008)", "MAT (2014/2015)", "SAB", "MIS", "APA", "CEK", "CHA", "IND", "HAR", "WAS", "SCA")


# Plot the PCA and save to file
temp_pca_plot <- ggplot(data=pca.dat) +
  geom_point(aes(x=Axis1, y=Axis2, colour=pop), size=10, alpha=0.8) +
  scale_colour_manual(values = col,
                        labels = labels) +
  rd_theme +
  theme(legend.title=element_blank(), legend.text=element_text(face="bold", size=14)) +
  labs(x = paste("", "Axis 1: ", eig_percent[1], "%"), y = paste("", "Axis 2: ", eig_percent[2], "%"))


ggsave(filename = "out/fig/figure_S1.png", plot = temp_pca_plot, width = 36, height = 24, units = "cm") 


