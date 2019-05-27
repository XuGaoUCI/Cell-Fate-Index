# This script is to reproduce Figure 5 and Figure 5S.
library(data.table)
library(openxlsx)
library(dplyr)
library(ggplot2)
library(stringi)
library(splines)
# load in data and preprocessing
human.signal <- read.xlsx('/path/to/human/signals.xlsx', sheet = 1)
rat.signal <- read.xlsx('/path/to/mice/signals.xlsx', sheet = 3)
cluster.human <- read.xlsx('/path/to/human/group.xlsx', sheet = 1)
cluster.rat <- read.xlsx('/path/to/mice/group.xlsx', sheet = 2)
roadmap <- read.xlsx('/path/to/roadmap.xlsx', sheet = 5)
colnames(roadmap) <- roadmap[1, ]
dic.human <- roadmap[1, 2:11]
dic.rat <- roadmap[1, 12:ncol(roadmap)]
# load in the library functions
source('/.../LibraryFunctions.R')
# generate the plots
for (name in dic.human) {
  GeneAnalysis(name, human = TRUE, scale = FALSE, output = FALSE, decreasing = FALSE)
}
# generate the plots
for (name in dic.human) {
  GeneAnalysis(name, human = TRUE, scale = TRUE, output = FALSE, decreasing = TRUE)
}
# generate the plots
for (name in dic.rat) {
  GeneAnalysis(name, human = FALSE, scale = FALSE, output = FALSE, decreasing = FALSE)
}
# generate the plots
for (name in dic.rat) {
  GeneAnalysis(name, human = FALSE, scale = TRUE, output = FALSE, decreasing = TRUE)
}
# overlaid human and mice 4 lists
over.list <- c("Contractility", "Cardiac structural protein", "Cell junction",
               "Ion channels")
human.cell <- GeneAnalysis(over.list, human = TRUE, scale = FALSE, output = TRUE, decreasing = FALSE)
rat.cell <- GeneAnalysis(over.list, human = FALSE, scale = FALSE, output = TRUE, decreasing = FALSE)
cell.name <- c(names(sort(human.cell)), names(sort(rat.cell)))
signal.combo <- c(sort(human.cell), sort(rat.cell))
species <- c(rep("Human", length(human.cell)), rep("Mice", length(rat.cell)))
cell.fraction <- c(1:length(human.cell)/length(human.cell), 1:length(rat.cell)/length(rat.cell))
over.df <- data.frame(cellfraction = cell.fraction , cellfateindex = signal.combo, 
                      SampleID = cell.name, species = species, stringsAsFactors = FALSE)
over.df <- over.df %>% inner_join(rbind(cluster.human[, 1:2], cluster.rat[, 1:2]), by = "SampleID")
ggplot(over.df, aes(x = cellfraction, y = cellfateindex, shape = species, group = species, color = GroupID)) + geom_point(size=1.5) + ylim(0, 1) + geom_smooth(method="auto", se = TRUE, level = 0.999) + 
  ggtitle("Human - Mice overlaid") + labs(x = "cell fraction", y = "cell fate index") +
  scale_color_manual(values=c("#619CFF", "#F8766D", "#C77CFF", "#7CAE00","#619CFF", "#F8766D", "#C77CFF")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        aspect.ratio = 1.5, text = element_text(size=15)) 
# overlaid plot speed and p value 
over.list <- c("Contractility", "Cardiac structural protein", "Cell junction",
               "Ion channels")
cell.matrix.rat <- c()
for (name in over.list) {
  cell.matrix.rat <- rbind(cell.matrix.rat, 
                           GeneAnalysis(name, human = FALSE, scale = FALSE, output = TRUE, decreasing = FALSE))
}

cell.matrix.human <- c()
for (name in over.list) {
  cell.matrix.human <- rbind(cell.matrix.human, 
                             GeneAnalysis(name, human = TRUE, scale = FALSE, output = TRUE, decreasing = FALSE))
}

cell.name <- c(as.vector(apply(cell.matrix.human, 1, function(x) names(sort(x)))), 
               as.vector(apply(cell.matrix.rat, 1, function(x) names(sort(x)))))
signal.combo <- c(as.vector(apply(cell.matrix.human, 1, function(x) (sort(x)))), 
                  as.vector(apply(cell.matrix.rat, 1, function(x) (sort(x)))))

species <- c(rep("Human", length(as.vector(cell.matrix.human))), rep("Mice", length(as.vector(cell.matrix.rat))))
cell.fraction <- c(rep(1:ncol(cell.matrix.human)/ncol(cell.matrix.human), 4), rep(1:ncol(cell.matrix.rat)/ncol(cell.matrix.rat), 4))
over.df <- data.frame(cellfraction = cell.fraction , cellfateindex = signal.combo, 
                      SampleID = cell.name, species = species, stringsAsFactors = FALSE)
over.df <- over.df %>% inner_join(rbind(cluster.human[, 1:2], cluster.rat[, 1:2]), by = "SampleID")
# scatterplot
aa <- qqplot(over.df[over.df$species == "Human", "cellfateindex"], 
             over.df[over.df$species == "Mice", "cellfateindex"], 
             xlab = "Human relative speed", ylab = "Mice relative speed")
abline(0,1)
# smoothing spline with separate absolute speed
fit.smooth.mice <- smooth.spline(1:length(aa$y)/length(aa$y),aa$y)
speed.mice <- CalculateDeri(fit.smooth.mice$x, fit.smooth.mice$y)
fit.smooth.human <- smooth.spline(1:length(aa$x)/length(aa$x),aa$x)
speed.human <- CalculateDeri(fit.smooth.human$x, fit.smooth.human$y)
over.speed <- data.frame(speed = c(speed.human, speed.mice), species = c(rep("human", length(speed.human)), rep("mice", length(speed.mice))), time = c(1:length(speed.human)/length(speed.human), 1:length(speed.human)/length(speed.human)))
ggplot(over.speed, aes(x = time, y = speed, color = species)) + geom_point(size=.1) + ylim(0, 1) + geom_smooth(method="auto", se = TRUE, level = 0.999) + 
  ggtitle("Human - Mice speed overlaid") + labs(x = "relative time", y = "speed") +
  scale_color_manual(values=c("#619CFF", "#F8766D")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        aspect.ratio = 1.5, text = element_text(size=15)) 

scatter.mice <- data.frame(time = fit.smooth.mice$x, speed = fit.smooth.mice$y)
scatter.human <- data.frame(time = fit.smooth.human$x, speed = fit.smooth.human$y)
ggplot(scatter.mice, aes(x = time, y = speed)) + geom_point(size=2.5, shape=c(2)) + ylim(0, 1) + geom_line(size = 0.8, color = 'red') + 
  ggtitle("Mice") + labs(x = "relative time", y = "cell fate index") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        aspect.ratio = 1.5, text = element_text(size=15)) 

ggplot(scatter.human, aes(x = time, y = speed)) + geom_point(size=2.5, shape=c(2)) + ylim(0, 1)  + geom_line(size = 0.8, color = 'red') + 
  ggtitle("Human") + labs(x = "relative time", y = "cell fate index") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        aspect.ratio = 1.5, text = element_text(size=15)) 
