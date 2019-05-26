# this script is to generate Figure 7 and Figure 7S.
library(data.table)
library(openxlsx)
library(dplyr)
library(ggplot2)
library(stringi)
# load in data and preprocessing
human.signal <- read.xlsx('/path/to/human/signals.xlsx', 
                          sheet = 1)
cluster.human <- read.xlsx('/path/to/mice/signals.xlsx', 
                           sheet = 2)
roadmap <- read.xlsx('path/to/roadmap.xlsx', 
                     sheet = 3)
source('/.../LibraryFunctions.R')
# generate the plots
GeneAnalysisDays(list.name = 'Cardiac structural protein', scale = FALSE, output = FALSE, decreasing = FALSE,
             is_filter = TRUE, 
             cr_slicer = c(1,2), cr_hccluster = c("ctrl","alternative", "reprogramming"),
             exclude = FALSE, shape.num = 16, genelist = 'Cardiac structural protein')
GeneAnalysisDays(list.name = 'Cardiac structural protein', scale = FALSE, output = FALSE, decreasing = FALSE,
             is_filter = TRUE, 
             cr_slicer = c(1,3), cr_hccluster = c("ctrl","alternative", "refractory"),
             exclude = FALSE, shape.num = 15, genelist = 'Cardiac structural protein')
GeneAnalysisDays(list.name = 'Contractility', scale = FALSE, output = FALSE, decreasing = FALSE,
             is_filter = TRUE, 
             cr_slicer = c(1,2,3), cr_hccluster = c( "ctrl", "reprogramming"),
             exclude = FALSE, shape.num = 16, genelist = 'Contractility')
GeneAnalysisDays(list.name = 'Contractility', scale = FALSE, output = FALSE, decreasing = FALSE,
             is_filter = TRUE, 
             cr_slicer = c(1,2,3), cr_hccluster = c( "ctrl","alternative"),
             exclude = FALSE, shape.num = 15, genelist = 'Contractility')
GeneAnalysisDays(list.name = 'Cell junction', scale = FALSE, output = FALSE, decreasing = FALSE,
             is_filter = TRUE, 
             cr_slicer = c(1,2), cr_hccluster = c("ctrl","alternative", "reprogramming"),
             exclude = FALSE, shape.num = 16, genelist = 'Cell junction')
GeneAnalysisDays(list.name = 'Cell junction', scale = FALSE, output = FALSE, decreasing = FALSE,
             is_filter = TRUE, 
             cr_slicer = c(1,3), cr_hccluster = c("ctrl","alternative", "refractory"),
             exclude = FALSE, shape.num = 15, genelist = 'Cell junction')
GeneAnalysisDays(list.name = 'ECM', scale = TRUE, output = FALSE, decreasing = TRUE,
             is_filter = TRUE, 
             cr_slicer = c(1,2,3), cr_hccluster = c("ctrl", "reprogramming"),
             exclude = FALSE, shape.num = 16, genelist = 'ECM', hcf = TRUE, y.lower= .4, y.upper = 1)
GeneAnalysisDays(list.name = 'ECM', scale = TRUE, output = FALSE, decreasing = TRUE,
             is_filter = TRUE, 
             cr_slicer = c(1,2,3), cr_hccluster = c("ctrl","alternative"),
             exclude = FALSE, shape.num = 15, genelist = 'ECM', y.lower= .4, y.upper = 1,hcf = TRUE)
GeneAnalysisDays(list.name = 'mRNA splicing', scale = TRUE, output = FALSE, decreasing = TRUE,
             is_filter = TRUE, 
             cr_slicer = c(1,2,3), cr_hccluster = c("ctrl", "reprogramming"),
             exclude = FALSE, shape.num = 16, genelist = 'mRNA splicing', hcf = TRUE, y.lower= .4, y.upper = 1)
GeneAnalysisDays(list.name = 'mRNA splicing', scale = TRUE, output = FALSE, decreasing = TRUE,
             is_filter = TRUE, 
             cr_slicer = c(1,2,3), cr_hccluster = c("ctrl","alternative"),
             exclude = FALSE, shape.num = 15, genelist = 'mRNA splicing', y.lower= .4, y.upper = 1,hcf = TRUE)
GeneAnalysisDays(list.name = 'Ion channels', scale = FALSE, output = FALSE, decreasing = FALSE,
             is_filter = TRUE, 
             cr_slicer = c(1,2,3), cr_hccluster = c( "ctrl", "reprogramming"),
             exclude = FALSE, shape.num = 16, genelist = 'Ion channels')
GeneAnalysisDays(list.name = 'Ion channels', scale = FALSE, output = FALSE, decreasing = FALSE,
             is_filter = TRUE, 
             cr_slicer = c(1,2,3), cr_hccluster = c( "ctrl","alternative"),
             exclude = FALSE, shape.num = 15, genelist = 'Ion channels')
#Supp F7
GeneAnalysisDays(list.name ='Cardiac_genelist_combined', scale = FALSE, output = FALSE, decreasing = FALSE,
             is_filter = TRUE, 
             cr_slicer = c(1,2,3), cr_hccluster = c("ctrl", "reprogramming"),
             exclude = FALSE, shape.num = 16,  genelist ='Cardiac_genelist_combined', all = TRUE)
GeneAnalysisDays(list.name = 'Cardiac_genelist_combined', scale = FALSE, output = FALSE, decreasing = FALSE,
             is_filter = TRUE, 
             cr_slicer = c(1,2,3), cr_hccluster = c("ctrl","alternative"),
             exclude = FALSE, shape.num = 15,  genelist ='Cardiac_genelist_combined', all = TRUE)
GeneAnalysisDays(list.name ='Fibroblast_genelist_combined', scale = TRUE, output = FALSE, decreasing = TRUE,
             is_filter = TRUE, 
             cr_slicer = c(1,2,3), cr_hccluster = c("ctrl", "reprogramming"),
             exclude = FALSE, shape.num = 16,  genelist ='Fibroblast_genelist_combined', all = TRUE,
             y.lower= .4, y.upper = 1, hcf = TRUE)

GeneAnalysisDays(list.name = 'Fibroblast_genelist_combined', scale = TRUE, output = FALSE, decreasing = TRUE,
             is_filter = TRUE, 
             cr_slicer = c(1,2,3), cr_hccluster = c("ctrl","alternative"),
             exclude = FALSE, shape.num = 15,  genelist = 'Fibroblast_genelist_combined', all = TRUE,
             y.lower= .4, y.upper = 1, hcf = TRUE)

