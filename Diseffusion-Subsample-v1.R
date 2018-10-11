# Project: Diseffusion: Subsample cells.
# Author: Zinaida Good
# Contact: zinaida@stanford.edu
# Version: v1 (10 October 2017)

# This script is designed to subsample desired numbaer of cells.
# Outputs: Modified sampled FCS file.
# Contributions: Current R implementation by Zinaida Good.

# Clear objects from the workspace.
rm(list = ls())

############################### Libraries ###############################
#source("http://bioconductor.org/biocLite.R") # Connect to BioConductor
#biocLite("flowCore") # Update BioConductor library for flow data

library(flowCore)  # Core package for FCS data
library(RcppCNPy)  # C++ package for 68K PBMC data.

########################### Define Constants ############################
# Number of cells to sample from an FCS file.
kSampleMin <- 1000
kSampleMid <- 10000
kSampleMax <- 50000

########################### Define Functions ############################
# Define a function to write an FCS file (based on Rachel's write.FCS function on GutHub).
cytofCore.write.FCS <- function(x, filename, what = "numeric", channelDescriptions = NULL, 
                                referenceDescription = NULL, oldDescription = NULL) {
  if (is.matrix(x)) {
    # Don't write "Leading_Push" to FCS files    
    x <- x[, which(colnames(x) != "Leading_Push")]
    
    # Build metadata for FCS file
    pd <- c()  # 'params' phenoData
    dl <- list()  # 'description' list
    
    dl[["$DATATYPE"]] <- "F"
    
    if (!is.null(referenceDescription)) {
      if (!is.null(referenceDescription[["$DATE"]])){
        dl[["$DATE"]] <- referenceDescription[["$DATE"]]
      }
      if (!is.null(referenceDescription[["$BTIM"]])){
        dl[["$BTIM"]] <- referenceDescription[["$BTIM"]]
      }
      if (!is.null(referenceDescription[["$ETIM"]])){
        dl[["$ETIM"]] <- referenceDescription[["$ETIM"]]
      }
      if (!is.null(referenceDescription[["$CYT"]])){
        dl[["$CYT"]] <- referenceDescription[["$CYT"]]
      }
      if (!is.null(referenceDescription[["$CYTSN"]])){
        dl[["$CYTSN"]] <- referenceDescription[["$CYTSN"]]
      }      
    }
    
    for (c in 1:ncol(x)) {
      c_name <- colnames(x)[c]
      c_desc <- colnames(x)[c]
      if (!is.null(channelDescriptions)){
        if (!is.na(channelDescriptions[c])){
          c_desc <- channelDescriptions[c]  
        }
      }
      
      c_min <- floor(min(0, min(x[, c])))  # Hack to prevent flowCore from shifting range
      c_max <- ceiling(max(x[, c]))
      c_rng <- c_max - c_min + 1
      
      pl <- matrix(c(c_name, c_desc, c_rng, c_min, c_max), nrow = 1)
      colnames(pl) <- c("name", "desc", "range", "minRange", "maxRange")
      rownames(pl) <- paste("$P", c, sep = "") 
      pd <- rbind(pd, pl)
      
      if (!is.null(referenceDescription[[paste0("P", c, "DISPLAY")]])){
        dl[[paste0("P", c, "DISPLAY")]] <- referenceDescription[[paste0("P", c, "DISPLAY")]]
      } 
      
      if (!is.null(referenceDescription[[paste0("$P", c, "G")]])){
        dl[[paste0("$P", c, "G")]] <- referenceDescription[[paste0("$P", c, "G")]]
      } 
      
      if (!is.null(referenceDescription[[paste0("$P", c, "R")]])){
        dl[[paste0("$P", c, "R")]] <- referenceDescription[[paste0("$P", c, "R")]]
      } else {
        dl[[paste("$P", c, "R",sep = "")]] <- toString(c_rng); # Range
      }
      
      if (!is.null(referenceDescription[[paste0("$P", c, "B")]])){
        dl[[paste0("$P", c, "B")]] <- referenceDescription[[paste0("$P", c, "B")]]
      } else {
        dl[[paste("$P", c, "B", sep = "")]] <- "32";      # Number of bits
      }
      
      if (!is.null(referenceDescription[[paste0("$P", c, "E")]])){
        dl[[paste0("$P", c, "E")]] <- referenceDescription[[paste0("$P", c, "E")]]
      } else { 
        dl[[paste("$P", c, "E", sep = "")]] <- "0,0";      # Exponent
      }
      
      dl[[paste("$P", c, "N", sep = "")]] <- c_name;      # Name
      dl[[paste("$P", c, "S", sep = "")]] <- c_desc;      # Desc	
    }	
    
    
    if (!is.null(oldDescription)) {
      if (!is.null(oldDescription[["$CYT"]])){
        dl[["$CYT"]] <- oldDescription[["$CYT"]]
      }		  
      if (!is.null(oldDescription[["$DATE"]])){
        dl[["$DATE"]] <- oldDescription[["$DATE"]]
      }
      if (!is.null(oldDescription[["$BTIM"]])){
        dl[["$BTIM"]] <- oldDescription[["$BTIM"]]
      }
      if (!is.null(oldDescription[["$ETIM"]])){
        dl[["$ETIM"]] <- oldDescription[["$ETIM"]]
      }
    }
    
    x <- flowFrame(x, as(data.frame(pd), "AnnotatedDataFrame"), description = dl) 
  }
  write.FCS(x, filename, what)
}

# Given an FCS file name, define a function to read in the data, 
# sample the desired number of cells, and write out as a new FCS file.
sampleFile <- function(file.name, data.set, n.sample) {  
  cat(paste("Reading :", file.name, "\n", sep = " "))
  file.fcs <- read.FCS(file.name, transformation = FALSE, emptyValue = FALSE)
  fcs.expr <- exprs(file.fcs)
  
  ## Get ground truth for Panorama.
  if (data.set == "Panorama") {
    annotation <- read.csv("Population_assignments.csv")[, c(1:2)]
    annotation <- annotation[order(annotation$Index.in.File), ]
    annotation <- annotation$Population
    annotation.df <- data.frame(Population = levels(factor(annotation)), 
                                n = table(as.numeric(factor(annotation))))
    fcs.expr <- data.frame(fcs.expr, ComponentID = factor(as.numeric(annotation)))
    #fcs.expr <- fcs.expr[fcs.expr$ComponentID != 1, ]
  }
  
  ## Sample the desired number of cells.
  if (nrow(fcs.expr) < n.sample) {
    sample.ids <- seq(1, nrow(fcs.expr))
  } else {
    sample.ids <- sample(seq(1, nrow(fcs.expr)), size = n.sample, replace = FALSE)
  }
  fcs.expr <- fcs.expr[sample.ids, ]
  
  ## Export ground truth summary for Panorama.
  if (data.set == "Panorama") {
    annotation.df <- data.frame(annotation.df,
                                Sampled = table(fcs.expr$ComponentID))
    write.csv(annotation.df, paste0("Population_legend_", n.sample, ".csv"))
  }
  
  ## Write out a modified FCS file containing sampled data.
  fcs.name <- paste0(substr(file.name, start = 1, stop = nchar(file.name) - 4), 
                     "_Sampled_", length(sample.ids), ".fcs")
  channels <- pData(parameters(file.fcs))[, "desc"]
  if (data.set == "Panorama") fcs.expr$ComponentID <- as.numeric(fcs.expr$ComponentID)
  cat(paste("Writing", fcs.name, "\n"))
  setwd(paste0(homeDir, "/", data.set, "/source_files"))
  suppressWarnings(cytofCore.write.FCS(as.matrix(fcs.expr), fcs.name, 
                                       what = "numeric", 
                                       channelDescriptions = channels,
                                       referenceDescription = description(file.fcs)))
}

################################# Inputs ################################
# 1. Select directories.
homeDir <- "~/Dropbox/Diseffusion/FCS-Inputs"
data.sets <- c("Panorama", "Wanderlust", "68K_PBMC", "Synth_10_50", "Synth_bifurcation")

############################# Process Data ##############################
# 1. Sample 1000, 10000, or 50000 cells from each FCS file.
setwd(homeDir)
#data.set <- "Panorama"

for(data.set in data.sets) {
  cat(paste("*** Sampling:", data.set, "***", "\n"))
  sourceDir <- paste0(homeDir, "/", data.set, "/source_files/all_cells")
  #sourceDir <- paste0(homeDir, "/", data.set, "/source_files/sample_50k")
  setwd(sourceDir)
  file.name <- list.files(pattern = ".fcs")[1]
  sampleFile(file.name, data.set, n.sample = kSampleMin)
  setwd(sourceDir)
  sampleFile(file.name, data.set, n.sample = kSampleMid)
  setwd(sourceDir)
  sampleFile(file.name, data.set, n.sample = kSampleMax)
}

# 2. Process 68K matrix.
setwd("~/Dropbox/Diseffusion/RNAseq-68K-PBMCs/fromCaleb")
npy68K <- npyLoad("PCA.maxMt_0.15.counts_3in3.minLogV_0.1.p_50.npy")
colnames(npy68K) <- vapply(1:ncol(npy68K), function(x) paste0("PC", x), "")
npy68K <- as.data.frame(npy68K)
setwd(sourceDir)
suppressWarnings(cytofCore.write.FCS(as.matrix(npy68K),
                                     filename = "RNAseq-68K-PBMCs.fcs", 
                                     what = "numeric", 
                                     channelDescriptions = colnames(npy68K)))

setwd(sourceDir)
file.name = "RNAseq-68K-PBMCs.fcs"
sampleFile(file.name, n.sample = kSampleMin)
setwd(sourceDir)
sampleFile(file.name, n.sample = kSampleMax)



