# library functions
TransformFunction <- function(signals, human = TRUE, scale = FALSE, start, end) {
  # this function is to do linear transform to the raw signals so that the starting point is 0
  # and the ending point is 1
  # Args 
  #   signals: the original data
  #   human: a boolean parameter to indicate human or mice
  #   scale: a boolean parameter to indicate whether to use start value as 0
  #   start: the initial list
  #   end: the ending list
  # Return
  #   transform.signals: the signals after transformation.
  
  signals <- unlist(signals)
  start.value <- mean(signals[start])
  end.value <- mean(signals[end])
  # determine the coefficients of linear equations
  linear.coef <- matrix(c(start.value, 1, end.value, 1), 2, 2, 
                        byrow = TRUE) %>% solve %*% c(0, 1)
  if (scale == TRUE) {
    linear.coef <- matrix(c(start.value, 1, end.value, 1), 2, 2, 
                          byrow = TRUE) %>% solve %*% c(1, 0)
  }
  signal.to.transform <- signals[!names(signals) %in% c(start, end)]
  transform.signals <- linear.coef[1] * signal.to.transform + linear.coef[2]
  # truncate within the range (0, 1)
  transform.signals[transform.signals <= 0] <- 0
  transform.signals[transform.signals >=1] <- 1
  return(transform.signals)
}

PlotSignalsDays <- function(transform.signals, title, decreasing = FALSE, shape.num = 19, y.lower= 0, y.upper = .4, hcf = FALSE) {
  # this function is to plot CFI by days.
  # Args 
  #   transform.signals: the transformed signals
  #   title: the title of the plot
  #   decreasing: a boolean parameter to indicate whether to use descending order
  #   shape.num: a paremter to choose the symbol of the plot
  #   y.lower: the minimum of y value
  #   y.upper: the maximum of y value
  #   hcf: a boolean parameter to indicate wheter to sort using descending or ascending order
  # Return
  #   g: a plot showing the cell fate index with cell fraction.
  
  transform.df.sort <- data.frame (signals = transform.signals, 
                                   Cell_ID = names(transform.signals), 
                                   stringsAsFactors = FALSE) %>%
    inner_join(cluster.human[, 1:2], by = "Cell_ID") 
  transform.df.sort['Days'] <- strtrim(transform.df.sort$Treatment, 2)
  if (hcf == FALSE) {
    transform.df.sort <- transform.df.sort %>%  arrange(Days, signals)
  } else {
    transform.df.sort <- transform.df.sort %>%  arrange(Days, desc(signals))
  }
  color.map <-  c("gold", "forestgreen", "mediumpurple", "red", "magenta3","black", 
                  "Blue","coral","darkred")
  names(color.map) <- c("D9M", "D7M", "D5M","D7R", "D5R", "D3UN","D3M", "D9R", "D3R")
  transform.df.sort <- transform.df.sort %>% cbind(index = (1:nrow(transform.df.sort))/nrow(transform.df.sort))
  pos <- transform.df.sort[,c('index', 'Days')] %>% group_by(Days) %>% summarise(pos = median(index))
  g <- ggplot(transform.df.sort, aes(x = index, y = signals)) + geom_point(shape = shape.num, size = 3,aes(colour = Treatment)) + ylim(y.lower, y.upper)  + 
    labs(x = "cell fraction", y = paste(title,"\ncell fate index")) +
    scale_color_manual(breaks = unique(transform.df.sort$Treatment),
                       values = color.map[unique(as.character(transform.df.sort$Treatment))]) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          aspect.ratio = .5, text = element_text(size=15),
          axis.text.x = element_text(size=8),
          axis.title.x=element_blank()) + 
    scale_x_discrete(limits=unlist(pos$pos),labels=unlist(pos$Days), expand = c(0.1,0))+
    guides(fill=guide_legend(title="Treatment"))
  print(g)
}

PlotSignals <- function(transform.signals, title, human = TRUE, decreasing = FALSE) {
  # this function is to plot CFI by days.
  # Args 
  #   transform.signals: the transformed signals
  #   title: the title of the plot
  #   human: a boolean parameter to indicate human or mice
  #   decreasing: a boolean parameter to indicate whether to use descending order
   # Return
  #   g: a plot showing CFI by cell fraction
  
  if (human == TRUE) {
    transform.df.sort <- data.frame (signals = transform.signals, 
                                     SampleID = names(transform.signals), 
                                     stringsAsFactors = FALSE) %>%
      inner_join(cluster.human[, 1:2], by = "SampleID") %>% 
      arrange(signals) 
    if (decreasing == TRUE) {
      transform.df.sort <- data.frame (signals = transform.signals, 
                                       SampleID = names(transform.signals), 
                                       stringsAsFactors = FALSE) %>%
        inner_join(cluster.human[, 1:2], by = "SampleID") %>% 
        arrange(desc(signals))
    }
    transform.df.sort <- transform.df.sort %>% cbind(index = (1:nrow(transform.df.sort))/nrow(transform.df.sort))
    g <- ggplot(transform.df.sort, aes(x = index, y = signals)) + geom_point(size = 3) + ylim(0, 1) + 
      aes(colour = GroupID) + ggtitle(title) + labs(x = "cell fraction", y = "cell fate index") +
      scale_color_manual(values=c("#619CFF", "#F8766D", "#C77CFF")) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            aspect.ratio = 1.5, text = element_text(size=15)) 
    print(g)
  } else {
    transform.df.sort <- data.frame (signals = transform.signals, 
                                     SampleID = names(transform.signals), 
                                     stringsAsFactors = FALSE) %>%
      inner_join(cluster.rat[, 1:2], by = "SampleID") %>% 
      arrange(signals) 
    if (decreasing == TRUE) {
      transform.df.sort <- data.frame (signals = transform.signals, 
                                       SampleID = names(transform.signals), 
                                       stringsAsFactors = FALSE) %>%
        inner_join(cluster.rat[, 1:2], by = "SampleID") %>% 
        arrange(desc(signals)) 
    }
    transform.df.sort <- transform.df.sort %>% cbind(index = (1:nrow(transform.df.sort))/nrow(transform.df.sort))
    g <- ggplot(transform.df.sort, aes(x = index, y = signals)) + geom_point(size = 3)  + ylim(0, 1) + 
      aes(colour = GroupID) + ggtitle(title) + labs(x = "cell fraction", y = "cell fate index") +
      scale_color_manual(values=c("#619CFF", "#F8766D", "#C77CFF", "#7CAE00")) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            aspect.ratio = 1.5, text = element_text(size=15))
    print(g)
  }
}

CalculateDeri <- function(domain, response) {
  # this function is to calculate the speed (derivatives) of CFI.
  # Args 
  #   domain: the values of x axis
  #   response: the values of y axis
  # Return
  #   deri: a list of speed of CFI
  
  length.curve <- length(domain)
  deri <- c()
  for (pos in 1:(length(domain) - 1)) {
    deri <- c(deri, (response[pos + 1] - response[pos]) / (domain[pos + 1] - domain[pos]))
  }
  return(deri)
}

GeneAnalysisDays <- function(list.name, scale = FALSE, output = FALSE,
                             decreasing = FALSE, is_filter = FALSE, 
                             cr_slicer = c(1,2,3), cr_hccluster = c("ctrl", "alternative", "reprogramming"),
                             exclude = FALSE, shape.num = 19, genelist, all = FALSE,  y.lower= 0, y.upper = .4, hcf = FALSE) {
  # the function is used to do the analyze and generate plots of CFI by days with different settings
  # Args 
  #   list.name: a gene list
  #   scale: a boolean parameter to indicate whether to use start value as 0
  #   output: a boolean paraeter to indicate whether to write the output
  #   decreasing: a boolean parameter to indicate whether to use descending order
  #   is_filter: a boolean parameter to show wheter to filter out some gene lists
  #   cr_slicer: a list of cluster of list
  #   cr_hccluster: a list of name of group of cells
  #   exclude: a boolean parameter to indicate whether to exclude some cells
  #   all:  a boolean parameter to indicate whether to include all gene lists
  #   shape.num: a paremter to choose the symbol of the plot
  #   y.lower: the minimum of y value
  #   y.upper: the maximum of y value
  #   hcf: a boolean parameter to indicate wheter to sort using descending or ascending order
  # Return
  #   transform.signals: transformed signals to write out and CFI plots
  
  roadmap <- roadmap[, 1:11]
  colnames(roadmap) <- roadmap[1,]
  roadmap <- roadmap[-c(1,2), -1]
  if (all) {
    if (hcf) {
      list.chosen <- roadmap[, 1:6] 
    } else {
      list.chosen <- roadmap[, 7:10] 
    }
  } else {
    list.chosen <- roadmap[,genelist]
  }
  list.chosen <- unlist(list.chosen)[!is.na(unlist(list.chosen))]
  filter.chosen <- unlist(cluster.human$Cell_ID)
  if (is_filter) {
    filter.chosen <- cluster.human %>% filter((cluster.human$SLICER_branch %in% cr_slicer) & 
                                                (cluster.human$HC_cluster %in% cr_hccluster)) %>% dplyr::select(Cell_ID)
    filter.chosen <- c('hCM1', 'hCM2', 'hCM3', unlist(filter.chosen))
    if (exclude) {
      exclude.chosen <- cluster.human$Cell_ID[grepl('D[0-9]R', cluster.human$Treatment)]
      filter.chosen <- filter.chosen[!filter.chosen %in% exclude.chosen]
    }
  }
  human.chosen <- human.signal %>% filter(human.signal[, 1] %in% list.chosen)
  row.names(human.chosen) <- human.chosen[, 1]
  if (is_filter | exclude) {
    human.chosen <- human.chosen[, unlist(filter.chosen)]
  } else {
    human.chosen <- human.chosen[,-1]
  }
  # log transform
  human.chosen <- log(human.chosen + 1) / log(2)
  start.index <- cluster.human %>% filter(cluster.human$Cell_ID %in% filter.chosen) %>% arrange((Pseudotime)) %>% select(Cell_ID) %>% head(10) %>% unlist
  end.index <- cluster.human$Cell_ID[1:3]
  print(start.index)
  if (decreasing == FALSE) {
    wgt <- human.chosen[, end.index] %>% rowMeans
    wgt <- wgt / sum(wgt)
  } else {
    wgt <- human.chosen[, start.index] %>% rowMeans
    wgt <- wgt / sum(wgt)
  }
  human.transform <- t(apply(human.chosen, 1, TransformFunction, scale = scale, start = start.index,
                             end = end.index))
  # weighted avg
  transform.signals <- wgt %*% human.transform %>% as.vector
  names(transform.signals) <- colnames(human.transform)
  PlotSignalsDays(transform.signals, title = list.name, decreasing = decreasing, shape.num = shape.num, 
                  y.lower= y.lower, y.upper = y.upper, hcf = hcf)
  if (output == TRUE) {
    return(transform.signals)
  }
}

GeneAnalysis <- function(list.name, human = TRUE, scale = FALSE, output = FALSE,
                         decreasing = FALSE) {
  # the function is used to do the analyze and generate plots of CFI with different settings
  # Args 
  #   list.name: a gene list
  #   scale: a boolean parameter to indicate whether to use start value as 0
  #   output: a boolean paraeter to indicate whether to write the output
  #   decreasing: a boolean parameter to indicate whether to use descending order
  #   human: a boolean parameter to indicate whether to use human or mice
  # Return
  #   transform.signals: transformed signals to write out and CFI plots
  
    if (human == TRUE) {
    roadmap <- roadmap[, 1:11]
    list.chosen <- roadmap[-c(1,2), list.name] 
    list.chosen <- unlist(list.chosen)
    human.chosen <- human.signal %>% filter(human.signal[, 1] %in% list.chosen)
    row.names(human.chosen) <- human.chosen[, 1]
    # log transform
    human.chosen <- log(human.chosen[, -1] + 1) / log(2)
    start.index <- cluster.human$SampleID[grepl("Start1|Start2", cluster.human[, 'X3'])]
    end.index <- cluster.human$SampleID[grepl("target cell", cluster.human[, 'X3'])]
    if (decreasing == FALSE) {
      wgt <- human.chosen[, end.index] %>% rowMeans
      wgt <- wgt / sum(wgt)
    } else {
      wgt <- human.chosen[, start.index] %>% rowMeans
      wgt <- wgt / sum(wgt)
    }
    human.transform <- t(apply(human.chosen, 1, TransformFunction, scale = scale, start = start.index,
                               end = end.index))
    # weighted avg
    transform.signals <- wgt %*% human.transform %>% as.vector
    names(transform.signals) <- colnames(human.transform)
    PlotSignals(transform.signals, title = paste("Human", list.name, sep = "-"), decreasing = decreasing)
  } else {
    roadmap <- roadmap[, 12:22]
    list.chosen <- roadmap[-c(1,2), list.name] 
    list.chosen <- unlist(list.chosen)
    rat.chosen <- rat.signal %>% filter(toupper(rat.signal[, 1]) %in% list.chosen)
    row.names(rat.chosen) <- rat.chosen[, 1]
    # log transform
    rat.chosen <- log(rat.chosen[, -1] + 1) / log(2)
    start.index <- cluster.rat$SampleID[grepl("Start1|Start2", cluster.rat[, 'X3'])]
    end.index <- cluster.rat$SampleID[grepl("target cell", cluster.rat[, 'X3'])]
    if (decreasing == FALSE) {
      wgt <- rat.chosen[, end.index] %>% rowMeans
      wgt <- wgt / sum(wgt)
    } else {
      wgt <- rat.chosen[, start.index] %>% rowMeans
      wgt <- wgt / sum(wgt)
    }
    rat.transform <- t(apply(rat.chosen, 1, TransformFunction, scale = scale, start = start.index,
                             end = end.index))
    # weighted avg
    transform.signals <- wgt %*% rat.transform %>% as.vector
    names(transform.signals) <- colnames(rat.transform)
    PlotSignals(transform.signals, title = paste("Mice", list.name, sep = "-"),
                human = FALSE, decreasing = decreasing)
  }
  if (output == TRUE) {
    return(transform.signals)
  }
}
