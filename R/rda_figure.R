library(adegenet)
library(ggplot2)
library(dplyr)
library(hierfstat)
library(vegan)
library(pegas)
library(grid)
library(gridExtra)


### Read input files ----------------------------

# Geographic coordinates (km - coastline distance)
coords = cbind(X=c(0, 290, 567, 1260, 1687, 1961, 2248, 2985, 3306, 3391, 3537))

# Make polynomials of spatial coordinates
space = poly(coords, degree=3)

# Give names to polynomials
colnames(space) <- c("X","X2","X3")

# Load the environmental data
tab <- read.table("data/raw/env_data.txt", header=TRUE, sep="\t", fill=TRUE, comment.char = "", stringsAsFactors=FALSE)
tab$Pop <- factor(tab$Pop, c("LLM", "MAT", "SAB", "MIS", "APA", "CEK", "CHA", "IND", "HAR", "WAS", "SCA"))
tab <- arrange(tab, Pop)
env <- tab[,3:ncol(tab)]
rownames(env) <- tab$Pop

# Remove non-continuous variables, geographic variables, and outliers
env <- select(env, -Estuary.Volume..m3., -Tide.Volume..m3., -Tide.Volume...Day..m3.d.1., -Estuary.Latitude, -Estuary.Longitude, -Catchment.Latitude, -Catchment.Longitude, -Tidal.Fresh.Blackwater, -Mixing.Zone.Blackwater, -Seawater.Blackwater, -Tides.day....)

# Change names with symbols to avoid problems with printing
env <- dplyr::rename(env, Oceanic.DIP = Oceanic.DIP..uM. , Oceanic.NO3 = Oceanic.NO3..uM.)

# Scale the variables to avoid problems with unequal variance
env <- as.data.frame(scale(env))


### RDA analysis - outlier loci ---------------------------

out_gen <- read.genepop("data/raw/final_outliers.gen", ncode = 3)

# Run a CA on the genind object
obj <- genind2genpop(out_gen)
ca1 <- dudi.coa(tab(obj),scannf=FALSE,nf=10)
barplot(ca1$eig,main="Correspondance Analysis eigenvalues",
        col=heat.colors(length(ca1$eig)))
s.label(ca1$li, sub="CA 1-2",csub=2)
out.data <- ca1$li[,1:10]
rownames(out.data) <- tab$Pop

# Perform forward selection of spatial variables - commented out to avoid rerunning

#ord.spa = rda(out.data ~ ., data.frame(space), scale= FALSE)
#stp.spa = ordistep(rda(out.data ~ 1, data.frame(space)), scope = formula(ord.spa), scale= FALSE, direction="forward", pstep = 1000)
#selected.spa = attributes(stp.spa$terms)$term.labels
selected_spa <- c("X", "X2")
space = space[,selected_spa]


# Perform forward selection of environmental variables - commented out to avoid rerunning

#ord.spa = rda(out.data ~ ., env, scale= FALSE)
#stp.spa = ordistep(rda(out.data ~ 1, env), scope = formula(ord.spa), scale= FALSE, direction="forward", pstep = 1000)
#selected.spa = attributes(stp.spa$terms)$term.labels
selected_env <- c("Oceanic.DIP", "Wind.Speed..m.sec.1.", "Ocean.Salinity.Min..psu.")
env = as.matrix(env[,selected_env])

# Perform the RDA
rda.out <- rda(out.data ~ space + env, scale=FALSE)
plot(rda.out)
RsquareAdj(rda.out)
vpart.out <- varpart(out.data, ~ space, ~ env)
plot(vpart.out)


# Create a biplot with genetic and environmental loadings
rda.out.eig <- rda.out$CCA$eig

rda.out.eig.perc <- (rda.out.eig / rda.out$CCA$tot.chi) * 100
rda.out.biplot <- as.data.frame(rda.out$CCA$biplot)
rda.out.biplot <- cbind(var = rownames(rda.out.biplot), rda.out.biplot)
rda.out.biplot <- filter(rda.out.biplot, var=="envOceanic.DIP" | var=="envWind.Speed..m.sec.1." | var=="envOcean.Salinity.Min..psu.")

rda.out.sites <- as.data.frame(rda.out$CCA$u)
rda.out.sites <- cbind(pop = rownames(rda.out.sites), rda.out.sites)

rda.out.axes <- rda.out$CCA$v

col <- c('#9e0142','#d53e4f','#f46d43','#fdae61','#a1d99b','#41ab5d','#00441b','#4eb3d3','#2b8cbe','#0868ac','#084081')
rda.out.sites$pop <- factor(rda.out.sites$pop, levels=rownames(rda.out.sites))

rda.out.sites <- mutate(rda.out.sites, region = case_when(pop %in% c("LLM", "MAT", "SAB", "MIS") ~ "NWG",
                                                        pop %in% c("APA", "CEK", "CHA") ~ "NEG",
                                                        pop %in% c("IND", "HAR", "WAS", "SCA") ~ "ATL")) %>%
  mutate(region = factor(region, levels = c("NWG", "NEG", "ATL")))

# Plot the biplot for the RDA
biplot <- ggplot(data=rda.out.sites) +
  geom_point(aes(x=RDA1, y=RDA2, colour=pop, shape = region), size=20, alpha=0.8) +
  scale_colour_manual(values = col,
                      labels = as.vector(rda.out.sites$pop),
                      name = "Locality") +
  scale_shape_manual(breaks = c("NWG", "NEG", "ATL"), values = c(16, 17, 18), name = "Region") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_segment(data=rda.out.biplot, aes(x=0, xend=RDA1, y=0, yend=RDA2), size=1.5, arrow = arrow(length = unit(0.025, "npc"))) +
  rd_theme +
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 18)) +
  labs(x = paste("", "Axis 1: ", sprintf("%.3f", rda.out.eig.perc[1]), "%"), y = paste("", "Axis 2: ", sprintf("%.3f", rda.out.eig.perc[2]), "%")) +
    xlim(-1, 1) +
    ylim(-0.65, 0.5) +
    coord_fixed(ratio=1) +
    annotate("text", x=-0.75, y=-0.65, label="Wind Speed", size=16) +
    annotate("text", x=-0.85, y=0.15, label="Oceanic DIP", size=16) +
    annotate("text", x=0.75, y=-0.65, label="Min Ocean Salinity", size=16)

# Get the relevant environmental variables

sig_vars_raw <- tab %>% select(Pop, Region, Oceanic.DIP..uM., Ocean.Salinity.Min..psu., Ocean.Salinity.Mean..psu., Ocean.Salinity.Max..psu., Wind.Speed..m.sec.1.)

sig_vars_raw$Region <- factor(sig_vars_raw$Region, levels = c("WGULF", "EGULF", "ATL"))

sal.plot <- ggplot(data=sig_vars_raw) +
  geom_point(aes(x=Pop, y=Ocean.Salinity.Mean..psu., colour = Region, shape = Region), size=10) +
    rd_theme +
    scale_colour_manual(name="Region",
                        breaks=c("WGULF", "EGULF", "ATL"),
                        labels=c("W. Gulf", "E. Gulf", "Atlantic"),
                        values=c("#9e0142", "#00441b", "#084081")) +
    scale_shape_manual(name="Region",
                     breaks=c("WGULF", "EGULF", "ATL"),
                     labels=c("NWG", "NEG", "ATL"),
                     values=c(16, 17, 18)) +
    ylab("Minimum Ocean Salinity") +
    xlab("Locality") +
    theme(axis.title.y = element_text(vjust = 0.5, angle = 90),
          axis.text.x = element_text(angle = 90),
          legend.position="none")

dip.plot <- ggplot(data=sig_vars_raw) +
  geom_point(aes(x=Pop, y=Oceanic.DIP..uM., colour = Region, shape = Region), size=10) +
  rd_theme +
  scale_colour_manual(name="Region",
                      breaks=c("WGULF", "EGULF", "ATL"),
                      labels=c("W. Gulf", "E. Gulf", "Atlantic"),
                      values=c("#9e0142", "#00441b", "#084081")) +
  scale_shape_manual(name="Region",
                     breaks=c("WGULF", "EGULF", "ATL"),
                     labels=c("NWG", "NEG", "ATL"),
                     values=c(16, 17, 18)) +
  ylab("Oceanic DIP (uM)") +
  xlab("Locality") +
  theme(axis.title.y = element_text(vjust = 0.5, angle = 90),
        axis.text.x = element_text(angle = 90),
        legend.position="none")

ws.plot <- ggplot(data=sig_vars_raw) +
  geom_point(aes(x=Pop, y=Wind.Speed..m.sec.1., colour = Region, shape = Region), size=10) +
  rd_theme +
  scale_colour_manual(name="Region",
                      breaks=c("WGULF", "EGULF", "ATL"),
                      labels=c("W. Gulf", "E. Gulf", "Atlantic"),
                      values=c("#9e0142", "#00441b", "#084081")) +
  scale_shape_manual(name="Region",
                     breaks=c("WGULF", "EGULF", "ATL"),
                     labels=c("NWG", "NEG", "ATL"),
                     values=c(16, 17, 18)) +
  ylab("Average Wind Speed (m/s)") +
  xlab("Locality") +
  theme(legend.position="none",
        axis.title.y = element_text(vjust = 0.5, angle = 90),
        axis.text.x = element_text(angle = 90))

legend <- get_legend(
    ggplot(data=sig_vars_raw) +
    geom_point(aes(x=Pop, y=Wind.Speed..m.sec.1., colour = Region, shape = Region), size=20) +
    rd_theme +
    scale_colour_manual(#name="Region",
                        breaks=c("WGULF", "EGULF", "ATL"),
                        labels=c("NWG", "NEG", "ATL"),
                        values=c("#9e0142", "#00441b", "#084081")) +
    scale_shape_manual(#name="Region",
                          breaks=c("WGULF", "EGULF", "ATL"),
                          labels=c("NWG", "NEG", "ATL"),
                          values=c(16, 17, 18)) +
    ylab("Average Wind Speed (m/s)") +
    xlab("Locality") +
    theme(axis.title.y = element_text(vjust = 0.5, angle = 90),
          axis.text.x = element_text(angle = 90),
          legend.key = element_blank(),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 26))
)


# Plot the figure

rda_plot <- grid.arrange(biplot, sal.plot, dip.plot, ws.plot, legend,
             layout_matrix = rbind(c(1,1,1, 1), c(2,3,4,5)),
             ncol=4, widths = c(8, 8, 8, 2), heights = c(16, 8))

ggsave(filename = "out/fig/figure_4.png", plot = rda_plot, width = 54, height = 48, units = "cm") 
