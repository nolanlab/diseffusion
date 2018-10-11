# Project: Diseffusion
# Author: Zinaida Good
# Contact: zinaida@stanford.edu
# Version: v5 (28 February 2018)

# Clear objects from the workspace.
rm(list = ls())

############################### Libraries ###############################
#source("http://bioconductor.org/biocLite.R")  # Connect to BioConductor
#biocLite("BiocUpgrade")  # Upgrate BioConductor installer
#biocLite("flowCore")  # Update BioConductor library for flowCore
#biocLite("FlowSOM")  # Update BioConductor library for FlowSOM

library(flowCore)  # Core package for FCS data
library(reshape2)  # Allows data format conversion with "melt"
#library(proxy)  # Similarity measures for 'dist()' including "cosine"
#library(MASS)  # Non-classical MDS (isoMDS)
#library(vegan)  # Isomap
#library(FlowSOM)  # Self-organizing maps
#browseVignettes("FlowSOM")
library(ggplot2)  # Good graphics
library(destiny)  # Diffusion maps
#browseVignettes("destiny")
library(RColorBrewer)  # Color pallets from RColorBrewer
#display.brewer.all()  # Can check available colors
library(vegan)  # Species diversity with "diversity"

########################### Define Constants ############################
kCells <- 1000  # Number of cells to assess in embedded space.
kMaxPlot <- 50000  # Max number of cells for plotting.
kScale <- FALSE  # NOT IMPLEMENTED # Scale data to mean = 0 and sd = 1 before embedding.
kDist <- "Euclidean"  # NOT IMPLEMENTED # Type of distance metric

########################### Define Functions ############################

################################ Inputs #################################
# 1. Select directories.
homeDir <- "/Users/zinaida/Documents/Stanford/Labs/Garry_P_Nolan/Diseffusion-Z"  # getwd()
srcDir <- paste0(homeDir, "/FCS-Inputs")
outDirName <- "Zina-Outputs"
setwd(homeDir)
#dir.create(outDirName)
outDir <- paste0(homeDir, "/", outDirName)

# 2. Select data sets and sample sets.
data.sets <- c("Synth_10_50", "Panorama", "68K_PBMC", "Synth_bifurcation", "Wanderlust")
sample.sets <- c("sample_1k", "sample_10k", "sample_50k")

############################# Process Data ##############################
#data.sets <- c("Synth_10_50", "Synth_bifurcation")
#sample.sets <- c("sample_1k", "sample_10k")
# data.set <- "Synth_bifurcation"
data.set <- "Synth_10_50"
sample.set <- "sample_1k"

# Process each dataset.
for(data.set in data.sets) {
  cat(paste("*** Processing:", data.set, "***", "\n"))
  
  # Set path to the directory containing source FCS file with viSNE and 
  # FDL, then process data.
  for(sample.set in sample.sets) {
    if (data.set == "Wanderlust") {
      if (sample.set == "sample_10k") sample.set <- "all_cells"
      if (sample.set == "sample_50k") break
    }
    cat(paste("***", data.set, "- Sample Set:", sample.set, "***", "\n"))
    
    sourceDir <- paste0(srcDir, "/", data.set, "/embeddings/", sample.set)
    setwd(outDir)
    dir.create(data.set)
    outputDir <- paste0(outDir, "/", data.set)
    setwd(outputDir)
    dir.create(sample.set)
    outputDir <- paste0(outDir, "/", data.set, "/", sample.set)
    
    source(paste0(homeDir, "/Diseffusion-Source-Process-File.R"))
  }
}
