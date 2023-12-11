set.seed(3112)
library(ggplot2)
library(tidyverse)

data <- read.csv("Nb_cells_per_cluster_Adjusted_Proportional.csv", header=T, row.names=1)

data_percentage <- apply(data, 2, function(x){x*100/sum(x,na.rm=T)})

pdf(file = "images/Nb_cells_per_cluster.pdf", width = 25, height = 4.5, pointsize=16.5)  # Adjust dimensions as needed

barplot(data_percentage,
        border = "black", 
        legend.text = TRUE, 
        args.legend = list(
          x = ncol(data_percentage) + 5.8,
          y = max(colSums(data_percentage)),
          bty = "n"
        ),
        ylab = "Percentage", 
        col = c("dodgerblue3", "tomato1", "palegreen4", "gold1"), 
        axisnames = TRUE,
        xlim = c(0, ncol(data_percentage) + 4.5),
        plot = TRUE)

dev.off()

