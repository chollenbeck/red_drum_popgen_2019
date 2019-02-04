library(ggplot2)
library(grid)


##### Plot themes and functions #####

# A custom ggplot theme for producing figures for the red drum popgen manuscript

rd_theme = theme(
  plot.title = element_text(size=16, face="bold", family="serif", vjust=1),
  panel.background = element_rect(fill = 'white'),
  panel.grid.major.y = element_line(size = 0.25, color = "lightgrey"),
  panel.grid.minor.y = element_blank(),
  panel.grid.major.x = element_line(size = 0.25, color = "lightgrey"),
  panel.grid.minor.x = element_blank(),
  axis.line = element_line(size = 0.7, color = "black"),
  axis.ticks = element_line(size = 0.25, color = "black"),
  axis.text = element_text(color="black", size=20, vjust=0.5, family="serif", angle = 0),
  axis.title.x = element_text(color="black", size=20, vjust=-0.35, face="bold"),
  axis.title.y = element_text(color="black", size=20, vjust=1, face="bold"),
  legend.key = element_rect(fill = "white"),
  panel.spacing = unit(0.25, "lines"),
  strip.background = element_rect(fill="white")
)


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

get_legend <- function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}