# Project: Diseffusion - Source to Process One FCS File
# Author: Zinaida Good
# Contact: zinaida@stanford.edu
# Version: v1 (28 February 2018)

############################### Libraries ###############################
require(flowCore)  # Core package for FCS data
require(ggplot2)  # Good graphics
require(destiny)  # Diffusion maps
require(reshape2)  # Allows data format conversion with "melt"
require(RColorBrewer)  # Color pallets from RColorBrewer
require(vegan)  # Species diversity with "diversity"

########################### Define Constants ############################
# Constants are defined in main.

########################### Define Functions ############################
# Define a function for arcsinh transformation (based on definition from Nikolay).
asinhNik <- function(value) {
  # value <- value - 1
  # for(i in 1:length(value)) {
  #   if((value[i] < 0) | is.na(value[i])) value[i] <- rnorm(1, mean = 0, sd = 0.01)
  # }
  value <- value / 5  # Co-factor value
  value <- asinh(value) # value <- log(value + sqrt(value^2 + 1))  
  return(value)
}

# Given an FCS file, define a function to read in data from the channels of interest.
get.channels <- function(file.fcs, arcsinh, xshift = FALSE) {
  ## Select channels of interest.
  channel.names <- pData(parameters(file.fcs))$desc
  channel.names <- as.vector(channel.names[grep(pattern = "(v)", channel.names)])
  channel.names <- c(channel.names, "tSNE1", "tSNE2")
  if ("ComponentID" %in% pData(parameters(file.fcs))$desc) {
    channel.names <- c(channel.names, "ComponentID")
  }
  if ("wanderlust" %in% pData(parameters(file.fcs))$desc) {
    channel.names <- c(channel.names, "wanderlust")
  }
  if (xshift) {
    channel.names <- tail(pData(parameters(file.fcs))$desc, 1)
  }
  ids <- which(pData(parameters(file.fcs))$desc %in% channel.names)
  channel.metals <- as.vector(pData(parameters(file.fcs))$name[ids])
  channels <- as.vector(pData(parameters(file.fcs))$desc[ids])
  x <- file.fcs[, channel.metals] 
  
  ## Arcsinh transform all marker channels except Division.
  if (arcsinh) {
    arcsinh.metals <- colnames(x)[-which(colnames(x) %in% 
                                           c("tSNE1", "tSNE2", "ComponentID", "wanderlust"))]
    arcsinht <- transformList(from = arcsinh.metals, tfun = asinhNik)
    x <- arcsinht %on% x
  }
  
  ## Create a data.frame containing transformed marker channels.
  x.df <- data.frame(exprs(x))
  if (!xshift) {
    rename.channels <- channels[-which(channels %in% 
                                         c("tSNE1", "tSNE2", "ComponentID", "wanderlust"))]
    channels[1:length(rename.channels)] <- as.vector(vapply(rename.channels, function(x) 
      substr(x, start = 1, stop = nchar(x) - 4), ""))
    names(x.df) <- channels
  }
  
  return(x.df)
}

# Given an FCS file name, define a function for read in data from that file, 
# and then extract the channels of interest for each cell.
get.fcs.data <- function(fcs.name, arcsinh, xshift = FALSE) {
  ## Extract and transform marker channels from the FCS file.
  file.fcs <- read.FCS(fcs.name, transformation = FALSE, emptyValue = FALSE)
  cat(paste("Reading", fcs.name, "\n"))
  fcs.df <- get.channels(file.fcs, arcsinh = arcsinh, xshift = xshift)
  rm(file.fcs)
  
  ## Scale data to mean = 0 and sd = 1 if needed.
  if (kScale) {
    fcs.df <- scale(fcs.df)
  }
  
  return(fcs.df)
}

# Define a function to calculate diffusion map and DPT.
getDiffusion <- function(fcs.expr, dist.metric = kDist) {
  if (dist.metric == "Euclidean") dist.metric <- "euclidean"
  ## Sample cells if file is too large.
  ncells <- nrow(fcs.expr)
  
  ## Angerer, Philipp (18May17): Calculate a lower bound approximation using 
  ## the number of memory bytes the matrix will need: 50000 * k * 16 
  ## (num cells * num nearest neighbors * bytes in a double-precision number *
  ## overhead of sparse matrix storage).
  k <- find_dm_k(ncells - 1)  # Default for n > 10k is k = 100, k = 1000 otherwise.
  if (ncells > 20000) k <- 20  # R crashes above this.

  ## Create and save a diffusion map.
  dm <- DiffusionMap(fcs.expr, distance = dist.metric, verbose = TRUE, k = k)
  saveRDS(dm, "Diffusion-Map-Object.rds")
  
  ## Plot diffusion component eigenvalues.
  picname <- paste0("Diffusion-Eigenvalues-", ncells, "-", dist.metric, ".eps")
  postscript(picname, width = 3.8, height = 4)
  plot(eigenvalues(dm), ylim = 0:1, pch = 20,
       xlab = 'Diffusion component (DC)', ylab = 'Eigenvalue')
  dev.off()
  
  ## Calculate Diffusion Pseudo-Time (DPT): a pseudo-time metric based on the
  ## transition probability of a diffusion process.
  # dpt <- DPT(dm)
  # saveRDS(dpt, "DPT-Object.rds")
  
  ## Plot DPT.
  # picname <- paste0("DPT-", ncells, "-", dist.metric, ".eps")
  # postscript(picname, width = 6, height = 5)
  # plot(dpt, pch = 20)
  # dev.off()
  
  ## Export cell coordinates in diffusion component space (not normalized).
  diffusion.df <- data.frame(DC1 = dm$DC1, DC2 = dm$DC2, DC3 = dm$DC3, DC4 = dm$DC4,
                             DC5 = dm$DC5, DC6 = dm$DC6, DC7 = dm$DC7, DC8 = dm$DC8,
                             DC9 = dm$DC9, DC10 = dm$DC10)
  # diffusion.df <- data.frame(diffusion.df,
  #                            DPT1 = dpt$DPT1, DPT2 = dpt$DPT2, DPT3 = dpt$DPT3,
  #                            DPT4 = dpt$DPT4, DPT5 = dpt$DPT5)
  write.csv(diffusion.df, "Diffusion-Cell-Coordinates.csv")
}

# Define a function to calculate Euclidean distance from the selected kCells 
# to all other cells in the embedded or original space. Export as CSV file.
getDistMx <- function(cells.df, distName) {
  cat(paste("Calculating", kDist, "distance in", distName, "space", "\n"))
  cells.df <- as.matrix(cells.df)
  
  if (nrow(cells.df) > 20000) {
    dist.df <- matrix(data = 0, nrow = nrow(cells.df), ncol = kCells)
    for(j in 1:kCells) {
      dist.df[, j] <- vapply(1:nrow(cells.df), 
                             function(i) dist(cells.df[c(i, j), ], 
                                              method = ifelse(kDist == "Euclidean", 
                                                              "euclidean", kDist))[1], 0)
    }
  } else {
    dist.df <- as.matrix(dist(cells.df), method = kDist)
  }

  write.csv(dist.df, paste0(distName, "-", kDist, "-Dist-Matrix.csv"))
}

# Define a function to plot data in embedded space.
plotData <- function(distName, annotation, colorScale = "zina1", plot.info = NULL) {
  ## Organize data for plotting.
  if (distName == "tSNE") {
    dist.df <- as.data.frame(fcs.data[, c("tSNE1", "tSNE2")])
  } else {
    dist.df <- read.csv(paste0(distName, "-Cell-Coordinates.csv"))[, -1]
    dist.df <- dist.df[, 1:2]
    if (distName == "PCA") colnames(dist.df) <- c("PC1", "PC2")
  }
  dist.df$Annotation <- annotation
  
  ## Plot data.
  colorName <- ""
  if (data.set == "Panorama") colorName <- "Population ID"
  if (data.set == "Synth_10_50" || data.set == "Synth_bifurcation") {
    colorName <- "Component ID"
  }
  if (data.set == "Wanderlust") colorName <- "Wanderlust"
  if (colorScale == "diversity") {
    colorName <- "Diversity"
    ## Bring cells with highest values to front.
    dist.df <- dist.df[order(dist.df$Annotation), ]
  }
  picname <- paste0(distName, "-by-Annotation", "-", colorScale, "-2D.eps")
  graphics.off()
  if (distName == "PCA") {
    p <- ggplot(dist.df, aes(x = PC1, y = PC2, colour = Annotation))
  } else if (distName == "tSNE") {
    p <- ggplot(dist.df, aes(x = tSNE1, y = tSNE2, colour = Annotation))
  } else if (distName == "scFDL") {
    p <- ggplot(dist.df, aes(x = X, y = Y, colour = Annotation))
  } else if (distName == "DC1-DC2") {
    p <- ggplot(dist.df, aes(x = DC1, y = DC2, colour = Annotation))
  } else if (distName == "DPT1-DPT2") {
    p <- ggplot(dist.df, aes(x = DPT1, y = DPT2, colour = Annotation))
  } else if (distName == "MDS") {
    p <- ggplot(dist.df, aes(x = MDS1, y = MDS2, colour = Annotation))
  } else if (distName == "IsoMDS") {
    p <- ggplot(dist.df, aes(x = IsoMDS1, y = IsoMDS2, colour = Annotation))
  } else if (distName == "Isomap") {
    p <- ggplot(dist.df, aes(x = Isomap1, y = Isomap2, colour = Annotation))
  }
  p <- p + geom_point(size = 0.5, shape = 21)
  if (colorScale != "default") {
    if (is.factor(annotation)) {
      if (colorScale == "zina1") {
        mycolors <- c("#8B2500", "#008B45", "#41B6C4", "#68228B", "#A1DAB4", "#5D476B",
                      "#225EA8", "#8B4513", "#191970", "#DDA0DD", "#2F4F4F", "#8B636C",
                      brewer.pal(8, "Dark2"), brewer.pal(8, "Set1"))
        
      } else if (colorScale == "zina2") {
        mycolors <- c(brewer.pal(8, "Dark2"), brewer.pal(8, "Set1"), brewer.pal(4, "YlGnBu"))
      }
      mycolors <- mycolors[1:length(levels(annotation))]
      p <- p + scale_colour_manual(values = mycolors)
    } else {
      if (colorScale == "zina1") {
        p <- p + scale_colour_gradient2(low = "#990000", mid = "#999999", high = "#002299",
                                        midpoint = 500)
      } else if (colorScale == "zina2") {
        p <- p + scale_colour_gradient2(low = "darkred", mid = "grey50", high = "darkblue",
                                        midpoint = 500)
      } else if (colorScale == "diversity") {
        # p <- p + scale_colour_gradient2(low = "#666666", mid = "#883333", high = "#990000", 
        #                                 midpoint = 0.5)
        # p <- p + scale_colour_gradient2(low = "#000099", mid = "#999922", high = "#FFFF00", 
        #                                 midpoint = 0.5)
        p <- p + scale_colour_gradient2(low = "#002299", mid = "#999999", high = "#990000",
                                        midpoint = 0.5)
        if (!is.null(plot.info)) {
          picname <- paste0(plot.info[1], "/", distName, "-by-Annotation", "-", 
                            colorScale, "-", plot.info[2], ".eps")
        }
      }
    }
  }
  p <- p + labs(title = colorName)
  p <- p + theme_bw()
  ggsave(filename = picname, width = 5, height = 4) 
}

# Define a function to plot correlation or rank agreement between Eucledian distance 
# in original space and in embedded space for K nearest neighbors.
makeCorPlots <- function(orig.df, dist.df, distName, ids, k, nCells = 10, type = "distance") {
  cat(paste("Plotting: Original vs. Embedded for", distName, type, "K =", k, "\n"))
  
  for(i in 0:nCells) {
    if (i == 0) {
      plot.df <- data.frame(Original = c(as.matrix(orig.df[ids, ])), 
                            Embedded = c(as.matrix(dist.df[ids, ])))
    } else {
      plot.df <- data.frame(Original = orig.df[ids, i], Embedded = dist.df[ids, i])
    }
    if (type == "rank") plot.df$Original <- seq(1:k) / nrow(orig.df)
    
    if (nrow(plot.df) > kMaxPlot) {
      plot.df <- plot.df[sample(1:nrow(plot.df), size = kMaxPlot, replace = FALSE), ]
    }
    
    picname <- paste0(distName, "-Cor-to-Original-", type, "-K", k, "-",
                      ifelse(i == 0, paste0(kCells, "-AllCells"), 
                             paste0("Cell", i)), ".png")
    kAlpha <- 0.01 * kMaxPlot / nrow(plot.df)
    if (kAlpha > 1) kAlpha <- 1
    graphics.off()
    p <- ggplot(plot.df, aes(x = Original, y = Embedded))
    p <- p + geom_point(size = 0.5, alpha = kAlpha)
    p <- p + geom_density2d()
    p <- p + labs(title = distName)
    p <- p + theme_bw()
    p <- p + expand_limits(x = c(0, 1), y = c(0, 1))
    ggsave(filename = picname, width = 4, height = 4) 
  }
}

# Define a function to calculate correlation between Eucledian distance in original 
# space and in embedded space for K nearest neighbors, with an option to plot data.
getCor <- function(orig.df, distName, K = nrow(orig.df) - 1, makePlots = FALSE) {
  dir.create(paste0(distName, "-Plots"))
  
  if (length(which(K > nrow(orig.df) - 1)) > 0) {
    K[which(K > nrow(orig.df) - 1)] <- nrow(orig.df) - 1
  }
  K <- unique(K)
  mycor.df <- matrix(data = 0, nrow = kCells, ncol = length(K))
  colnames(mycor.df) <- K
  mycorPearson.df <- mycor.df
  median.rank <- mycor.df
  percent.rank <- mycor.df
  knn.density <- mycor.df
  
  ## Load embedded data.
  dist.df <- read.csv(paste0(distName, "-", kDist, "-Dist-Matrix.csv"))[, -1]
  if (distName == "DPT") dist.df <- data.frame(DPT = dist.df)
  
  ## Sort cells by nearest neignbors in the original space.
  order.df <- as.data.frame(apply(orig.df, MARGIN = 2, function(x) order(x)))
  for(i in 1:kCells) {
    dist.df[, i] <- dist.df[, i][order.df[, i]]
    orig.df[, i] <- orig.df[, i][order.df[, i]]
  }
  
  ## Remove reference to self.
  orig.df <- orig.df[-1, ]
  dist.df <- dist.df[-1, ]
  
  ## Define ranks for all neighbors in the embedded space (on line 0 - 1).
  rank.df <- as.data.frame(apply(dist.df, MARGIN = 2, function(x) rank(x)))
  
  ## Calculate Spearman correlation between original and embedded distances for each K.
  setwd(paste0(distName, "-Plots"))
  for (j in 1:length(K)) {
    k <- K[j]
    col.id <- which(K == k)
    ids <- 1:nrow(orig.df)
    if (k < nrow(orig.df)) ids <- 1:k
    
    ## Spearman orrelation.
    mycor <- vapply(1:kCells, function(i) cor(x = orig.df[ids, i], 
                                              y = dist.df[ids, i], method = "spearman"), 0)
    mycor.df[, col.id] <- mycor
    
    ## Pearson correlation.
    mycor <- vapply(1:kCells, function(i) cor(x = orig.df[ids, i], 
                                              y = dist.df[ids, i], method = "pearson"), 0)
    mycorPearson.df[, col.id] <- mycor
    
    ## Calculate median rank of cells in the embedded space (as a fraction of K).
    myrank <- vapply(1:kCells, function(i) median(rank.df[ids, i]), 0)
    median.rank[, col.id] <- myrank * (2 / k)
    
    ## Calculate percent of K nearest neighbors captured by K in the embedded space.
    myrank <- vapply(1:kCells, function(i) mean(rank.df[ids, i] <= k), 0)
    percent.rank[, col.id] <- myrank
    
    ## Calclate KNN density (0-1) in the original space: 1 / (1 + sum(dist)).
    mydensity <- vapply(1:kCells, function(i) 1 / (1 + sum(orig.df[ids, i])), 0)
    knn.density[, col.id] <- mydensity
    
    ## Plot results for K = k.
    if (makePlots) {
      makeCorPlots(orig.df, dist.df, distName, ids, k, nCells = 10, type = "distance")
      makeCorPlots(orig.df, rank.df, distName, ids, k, nCells = 10, type = "rank")
    }
  }
  setwd("../")
  
  ## Write correlation.
  write.csv(mycorPearson.df, paste0(distName, "-Distance-Pearson-Correlation.csv"))
  write.csv(mycor.df, paste0(distName, "-Distance-Spearman-Correlation.csv"))
  write.csv(median.rank, paste0(distName, "-Median-Cell-Rank-Fraction.csv"))
  write.csv(percent.rank, paste0(distName, "-Fraction-Cells-by-Rank.csv"))
  write.csv(knn.density, paste0(distName, "-KNN-Density.csv"))
}

# Define a function to read in correlation or rank results for a given distance metric.
readCor <- function(distName, type) {
  if (type == "median-rank") type <- "Median-Cell-Rank-Fraction"
  if (type == "rank-fraction") type <- "Fraction-Cells-by-Rank"
  if (type == "dist-pearson-cor") type <- "Distance-Pearson-Correlation"
  if (type == "dist-spearman-cor") type <- "Distance-Spearman-Correlation"
  if (type == "knn-density") type <- "KNN-Density"
  distName2 <- ifelse(distName == "Diffusion", "DC1-DC2", distName)
  cor.df <- read.csv(paste0(distName2, "-", type, ".csv"))[, -1]
  colnames(cor.df) <- vapply(colnames(cor.df), 
                             function(x) sub(pattern = "X", replacement = "K", x), "")
  #colnames(cor.df)[ncol(cor.df)] <- "Kmax"
  cor.df$Method <- rep(distName, nrow(cor.df))
  return(cor.df)
}

# Define a function to plot correlation results as a modified boxplot.
plotCor <- function(cor.df, type, violin = FALSE) {
  if (type == "rank-fraction") y.name <- "Fraction KNN captured"
  if (type == "median-rank") y.name <- "Median rank (* 2/K)"
  if (type == "dist-pearson-cor") y.name <- "Pearson Correlation: orig. vs. emb. dist."
  if (type == "dist-spearman-cor") y.name <- "Spearman Correlation: orig. vs. emb.dist."
  
  melt.df <- melt(cor.df, id.vars = "Method", variable.name = "K", 
                  value.name = "Correlation")
  melt.df$Method <- factor(melt.df$Method, 
                           levels = c("PCA", "Diffusion", "scFDL", "tSNE"))
  
  p <- ggplot(melt.df, aes(x = K, y = Correlation, fill = Method))
  if (violin) {
    dodge <- position_dodge(width = 0.9)
    p <- p + geom_violin(position = dodge, scale = "width", trim = TRUE)
    p <- p + geom_boxplot(outlier.shape = 1, outlier.size = 1,
                          position = dodge, width = 0.2)
  } else {
    p <- p + geom_boxplot(outlier.shape = 1, outlier.size = 1)
  }
  p <- p + scale_fill_brewer(type = "seq", palette = "YlGnBu")
  p <- p + expand_limits(y = c(0, 1))
  p <- p + labs(title = "", x = "Number of Nearest Neighbors", y = y.name)
  if (type == "median-rank") p <- p + scale_y_log10()
  p <- p + theme_bw()
  
  ## Export plot as EPS file.
  picname <- paste0("Original-vs-Embedded-Distances-", type, "-by-K-",
                    ifelse(violin, "Violin", "Boxplot"), ".eps")
  ggsave(picname, width = ifelse(violin, 5, 4), height = 4)
}

# Define a function to plot a given evaluation metric type as a function of 
# KNN density in the original space for a given k.
plotKnnDensity <- function(type) {
  ## Read in data for all K.
  metric.df <- read.csv(paste0("Embedded-vs-Original-", type, ".csv"), 
                        header = TRUE, stringsAsFactors = FALSE)[, -1]
  knn.df <- read.csv(paste0("Embedded-vs-Original-knn-density.csv"),
                     header = TRUE, stringsAsFactors = FALSE)[, -1]
  Ks <- colnames(knn.df)[-ncol(knn.df)]
  
  ## Plot data for each K.
  setwd(paste0(outputDir, "/KNN-Density-Plots"))
  for(k in Ks) {
    cat(paste("Plotting:", type, "over KNN density, K =", k, "\n"))
    
    ## Select data for a given K to plot.
    plot.df <- data.frame(x = knn.df[, k], y = metric.df[, k], Method = knn.df$Method)
    plot.df$x <- (plot.df$x - min(plot.df$x)) / (max(plot.df$x) - min(plot.df$x))
    plot.df$Method <- factor(plot.df$Method, levels = c("PCA", "Diffusion", "scFDL", "tSNE"))
    graphics.off()
    picname <- paste0("KNN-Density-", type, "-", k)
    y.name <- ""
    if (type == "rank-fraction") y.name <- "Fraction KNN captured"
    if (type == "median-rank") y.name <- "Median rank (* 2/K)"
    if (type == "dist-pearson-cor") y.name <- "Pearson Correlation: orig. vs. emb. dist."
    if (type == "dist-spearman-cor") y.name <- "Spearman Correlation: orig. vs. emb.dist."
    
    ## Plot data as a scatter plot.
    p <- ggplot(plot.df, aes(x = x, y = y, color = Method))
    p <- p + geom_point(size = 0.5, shape = 21)
    p <- p + labs(x = "KNN Density (%)", y = y.name, title = k)
    p <- p + theme_bw()
    p <- p + expand_limits(x = c(0, 1), y = c(0, 1))
    p <- p + scale_color_brewer(type = "seq", palette = "YlGnBu")
    ggsave(filename = paste0(picname, "-scatter.eps"), width = 4, height = 4) 
    
    ## Plot KNN density distribution as a histogram.
    png(filename = paste0(picname, "-hist.png"), width = 480, height = 480)
    hist(plot.df$x)
    dev.off()
    
    ## Plot data as a box plot.
    plot.df$x_buckets <- cut(plot.df$x,
                             breaks = quantile(plot.df$x, probs = seq(0, 1, by = 0.2), na.rm = TRUE),
                             include.lowest = TRUE)
    p <- ggplot(plot.df, aes(x = x_buckets, y = y, fill = Method))
    p <- p + geom_boxplot(outlier.shape = 1, outlier.size = 1)
    p <- p + labs(x = "KNN Density", y = y.name, title = k)
    p <- p + theme_bw()
    p <- p + expand_limits(y = c(0, 1))
    p <- p + scale_fill_brewer(type = "seq", palette = "YlGnBu")
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    ggsave(filename = paste0(picname, "-boxplot.eps"), width = 4, height = 4) 
  }
  
  setwd("../")
}

# Define a function to calculate distance between 2 cells scaled to global sigma.
# Used in trajectory continuity calculation.
getDist <- function(cell1, cell2, sigma = NULL) {
  dist.1.2 <- rbind(cell1, cell2)
  
  if(!is.null(sigma)) {
    for (j in 1:ncol(dist.1.2)) {
      dist.1.2[, j] <- dist.1.2[, j] / sigma[j]  # Normalize to global sigma.
    }
  }

  dist.1.2 <- as.matrix(dist(dist.1.2), method = kDist)[1, 2]
  
  return(dist.1.2)
}

# Define a function to plot trajectory continuity for a given distance metric.
plotTrajContinuity <- function(distName, annotation, colorScale = "zina1",
                               plot.dir = "Trajectory-Continuity", plotSD = FALSE) {
  ## Import cell coordinates in embedded space.
  if (distName == "tSNE") {
    coord.df <- as.data.frame(fcs.data[, c("tSNE1", "tSNE2")])
  } else if (distName == "Original") {
    coord.df <- as.data.frame(fcs.data[, -which(colnames(fcs.data) %in% c("tSNE1", "tSNE2"))])
  } else {
    coord.df <- read.csv(paste0(distName, "-Cell-Coordinates.csv"))[, -1]
    coord.df <- coord.df[, 1:2]
    if (distName == "PCA") colnames(coord.df) <- c("PC1", "PC2")
  }
  coord.df$Annotation <- annotation
  
  ## Calculate global SD for both embedded dimentions to normalize distance
  # to Z-score difference.
  sigma <- apply(coord.df[, -ncol(coord.df)], MARGIN = 2, sd)
  
  ## Calculate trajectory continuity (distance between neighboring bins).
  is <- sort(unique(annotation))
  dist.df <- data.frame(i1 = is[-length(is)],  # Bin 1
                        n1 = rep(-1, length(is) - 1),  # n cells in bin 1
                        i2 = is[-1],  # Bin 2
                        n2 = rep(-1, length(is) - 1),  # n cells in bin 2
                        Distance = rep(-1, length(is) - 1),  # Distance norm. to sigma
                        SD = rep(-1, length(is) - 1))  # SD of distance norm. to sigma
  
  ## Synth_bifurcation has 1000 components; 500 is the last prior to branching;
  ## 501, 503, etc. (odd) is branch 1; 502, 504, etc. (even) is branch 2.
  if (data.set == "Synth_bifurcation") {
    branches.df <- subset(dist.df, i1 > 500)
    i2s <- vapply(branches.df$i1, function(x) {
      if (x %% 2 == 0) {
        return(branches.df$i2[which(branches.df$i2 %% 2 == 0 & branches.df$i2 > x)[1]])
      } else {
        return(branches.df$i2[which(branches.df$i2 %% 2 == 1 & branches.df$i2 > x)[1]])
      }
    }, 0)
    dist.df$i2[dist.df$i1 > 500] <- i2s
  }
  dist.df <- dist.df[!is.na(dist.df$i2), ]
  
  for (i in 1:nrow(dist.df)) {
    ## Select bin 1 cells.
    i1 <- dist.df$i1[i]
    bin1 <- subset(coord.df, Annotation == i1)
    centroid1 <- apply(bin1[, -ncol(bin1)], MARGIN = 2, mean)
    
    ## Select bin 2 cells.
    i2 <- dist.df$i2[i]
    bin2 <- subset(coord.df, Annotation == i2)
    centroid2 <- apply(bin2[, -ncol(bin2)], MARGIN = 2, mean)
    
    ## Calculate distance.
    centroid.dist <- getDist(centroid1, centroid2, sigma)
    if (plotSD) {
      cells.dist <- matrix(data = 0, nrow = nrow(bin2), ncol = nrow(bin1))
      for(k in 1:nrow(bin2)) {
        cell2 <- bin2[k, ]
        cells.dist[k, ] <- apply(bin1, MARGIN = 1, function(cell1) getDist(cell1, cell2, sigma))
      }
      sd.dist <- sd(cells.dist)
    } else {
      sd.dist <- -1
    }
    dist.df[dist.df$i1 == i1, -1] <- c(nrow(bin1), i2, nrow(bin2), centroid.dist, sd.dist)
  }
  
  ## Filter for data on neighboring components only.
  keep.ids <- dist.df$i2 - dist.df$i1
  if (data.set == "Synth_bifurcation") {
    ## Correct for branching after 500.
    keep.ids[which(dist.df$i1 > 500)] <- keep.ids[which(dist.df$i1 > 500)] / 2
  }
  keep.ids <- which(keep.ids == 1)
  dist.df <- dist.df[keep.ids, ]
  
  ## Plot data as a dotplot.
  setwd(plot.dir)
  picname <- paste0(distName, "-Continuity")
  graphics.off()
  p <- ggplot(dist.df, aes(x = i1, y = Distance, colour = i1))
  p <- p + geom_point(size = 0.5, shape = 21)
  if (colorScale != "default") {
    if (colorScale == "zina1") {
      p <- p + scale_colour_gradient2(low = "#990000", mid = "#999999", high = "#002299",
                                      midpoint = 500)
    } else if (colorScale == "zina2") {
      p <- p + scale_colour_gradient2(low = "darkred", mid = "grey50", high = "darkblue",
                                      midpoint = 500)
    }
  }
  p <- p + labs(x = "Component ID", y = "Distance between components")
  p <- p + theme_bw()
  if (data.set == "Wanderlust") p <- p + expand_limits(y = c(0, 4))
  if (data.set == "Synth_bifurcation") p <- p + expand_limits(y = c(0, 0.6))
  ggsave(filename = paste0(picname, "-Dotplot.eps"), width = 5, height = 4)
  
  ## Plot data as a histogram.
  p <- ggplot(dist.df, aes(Distance))
  p <- p + geom_density(fill = "grey50", trim = FALSE)
  p <- p + labs(y = "Distance between components")
  p <- p + theme_bw()
  if (data.set == "Wanderlust") p <- p + expand_limits(x = c(0, 4))
  if (data.set == "Synth_bifurcation") p <- p + expand_limits(x = c(0, 0.6))
  ggsave(filename = paste0(picname, "-Histogram.eps"), width = 4, height = 4)
  
  ## Output distance data as a CSV file.
  write.csv(dist.df, paste0(distName, "-Continuity-Data.csv"))
  
  ## Calculate and output distance statistics.
  dist.stats <- data.frame(Median = median(dist.df$Distance), 
                           Mean = mean(dist.df$Distance), 
                           SD = sd(dist.df$Distance), 
                           as.data.frame(t(quantile(dist.df$Distance, 
                                                    probs = c(0.25, 0.5, 0.75, 0.9, 0.95, 0.99)))))
  if (data.set == "Wanderlust") {
    dist.stats <- data.frame(dist.stats,
                             P050 = mean(dist.df$Distance > 0.5),
                             P075 = mean(dist.df$Distance > 0.75),
                             P100 = mean(dist.df$Distance > 1.0))
  } else if (data.set == "Synth_bifurcation") {
    dist.stats <- data.frame(dist.stats,
                             P0025 = mean(dist.df$Distance > 0.025),
                             P0050 = mean(dist.df$Distance > 0.050),
                             P0100 = mean(dist.df$Distance > 0.100))
  }
  write.csv(dist.stats, paste0(distName, "-Continuity-Stats.csv"))
  
  ## Plot original vs. embedded continuity.
  if (distName != "Original" && plotSD) {
    orig.dist <- read.csv("Original-Continuity-Data.csv", 
                          header = TRUE, stringsAsFactors = FALSE)[, -1]
    plot.df <- data.frame(Original = cumsum(orig.dist$Distance),
                          Embedded = cumsum(dist.df$Distance),
                          SD = dist.df$SD,
                          ComponentID = dist.df$i1)
    plot.df$Original <- plot.df$Original / max(plot.df$Original)
    plot.df$SD <- plot.df$SD / max(plot.df$Embedded)
    plot.df$Embedded <- plot.df$Embedded / max(plot.df$Embedded)
    plot.df$SD[which(is.na(plot.df$SD))] <- 0
    plot.df$Lower <- plot.df$Embedded - plot.df$SD
    plot.df$Upper <- plot.df$Embedded + plot.df$SD
    
    p <- ggplot(plot.df, aes(x = Original, y = Embedded, colour = ComponentID))
    p <- p + geom_pointrange(aes(ymin = Lower, ymax = Upper), shape = 20, size = 0.1)
    # p <- p + geom_linerange(aes(ymin = Lower, ymax = Upper), size = 0.5)
    if (colorScale != "default") {
      if (colorScale == "zina1") {
        p <- p + scale_colour_gradient2(low = "#990000", mid = "#999999", high = "#002299",
                                        midpoint = 500)
      } else if (colorScale == "zina2") {
        p <- p + scale_colour_gradient2(low = "darkred", mid = "grey50", high = "darkblue",
                                        midpoint = 500)
      }
    }
    p <- p + labs(x = "Original", y = "Embedded")
    p <- p + theme_bw()
    ggsave(filename = paste0(picname, "-by-Component-Dotplot.eps"), width = 5, height = 4)
  }
  
  setwd("../")
}

# Define a function to calculate diversity index for each KNN in Ks in embedded space.
getPopPurity <- function(distName, annotation, Ks = Ks, plot.dir = "Population-Purity",
                         div.type = "Gini-Simpson") {
  ## Load embedded data.
  dist.df <- read.csv(paste0(distName, "-", kDist, "-Dist-Matrix.csv"))[, -1]
  
  ## Calculate Gini-Simpson index (1 - Simpson index), 1 - Shannon index, or
  ## Inverse Simpson (1 / Simpson index) for each KNN in embedded space.
  ## (Inverse Simpson index = Inf if all cells are from the same population).
  for (div.type in c("Gini-Simpson", "1-Shannon", "Inv-Simpson")) {
    if (div.type == "Gini-Simpson") div.index <- "simpson"
    if (div.type == "1-Shannon") div.index <- "shannon"
    if (div.type == "Inv-Simpson") div.index <- "invsimpson"
    gs.mx <- matrix(data = 0, nrow = nrow(dist.df), ncol = length(Ks))
    colnames(gs.mx) <- Ks
    for(i in 1:length(Ks)) {
      K = Ks[i]
      gs.index <- apply(dist.df, MARGIN = 2,
                        function(x) { diversity(as.numeric(annotation[order(x)][2:(K + 1)]),
                                                    index = div.index)} )
      gs.mx[, i] <- gs.index
    }
    if (div.type == "Gini-Simpson") gs.mx <- 1 - gs.mx
    
    ## Write diversity index data.
    write.csv(gs.mx, paste0(plot.dir, "/", distName, "-", div.type, "-Index.csv"))
  }
}

# Define a function to plot diversity index in embedded space after normalizing 
# data to 0 - 1 for each KNN in Ks among all embedding methods.
plotPopPurity <- function(distNames, annotation, K, plot.dir = "Population-Purity",
                          div.type = "Gini-Simpson", to.log = FALSE, scale.div = 1, clamp = 1) {
  setwd(plot.dir)
  
  ## Load diversity index data for a given K.
  gs.mx <- matrix(data = 0, nrow = nrow(fcs.data), ncol = length(distNames))
  colnames(gs.mx) <- distNames
  for(distName in distNames) {
    gs.mx[, distName] <- read.csv(paste0(distName, "-", div.type,
                                         "-Index.csv"))[, paste0("X", K)]
  }
  gs.mx <- gs.mx[, c(2, 4, 3, 1)]
  gs.min <- min(gs.mx)
  gs.max <- max(gs.mx)
  cat(paste(div.type, "K =", K, ": min =", min(gs.mx), "max =",  max(gs.mx), "\n"))
  
  ## Plot data as a histogram.
  for(distName in distNames) {
    plot.df <- as.data.frame(gs.mx[, distName], drop = FALSE)
    colnames(plot.df)[1] <- "Diversity"
    p <- ggplot(plot.df, aes(Diversity))
    p <- p + geom_density(fill = "grey50", trim = FALSE)
    p <- p + labs(y = paste0(div.type, " index"))
    p <- p + theme_bw()
    p <- p + expand_limits(x = c(gs.min, gs.max))
    ggsave(filename = paste0(div.type, "-Histogram-", distName, "-K", K, ".eps"), 
           width = 4, height = 4)
  }
  
  ## Plot data as a boxplot.
  plot.df <- as.data.frame(gs.mx)
  colnames(plot.df)[2:3] <- c("DM", "FDL")
  plot.df <- melt(plot.df, variable.name = "Method", value.name = "Diversity")
  p <- ggplot(plot.df, aes(x = Method, y = Diversity, fill = Method))
  p <- p + geom_boxplot(outlier.shape = 1, outlier.size = 1)
  p <- p + labs(x = "", y = paste0(div.type, " diversity index"))
  p <- p + theme_bw()
  p <- p + scale_fill_brewer(type = "seq", palette = "YlGnBu")
  ggsave(filename = paste0(div.type, "-K", K, "-Boxplot.eps"), width = 4, height = 4) 
  
  ## Output distance data as a CSV file.
  write.csv(gs.mx, paste0(div.type, "-Diversity-Data-K", K, ".csv"))
  
  ## Calculate and output distance statistics.
  dist.stats <- apply(gs.mx, MARGIN = 2, 
                      function(x) c(Mean = mean(x), SD = sd(x), Median = median(x), 
                                    Min = min(x), 
                                    quantile(x, probs = c(0.25, 0.5, 0.75, 0.9, 0.95, 0.99)),
                                    Max = max(x)))
  write.csv(dist.stats, paste0(div.type, "-Diversity-Stats-K", K, ".csv"))
  
  ## Normalize diversity index to 0-1 range; brighten low-value cells if needed.
  if (to.log) gs.mx <- log(gs.mx)
  gs.mx <- gs.mx - min(gs.mx)
  gs.mx <- gs.mx / max(gs.mx)
  if (div.type == "1-Shannon" || div.type == "Inv-Simpson") gs.mx <- 1 - gs.mx
  gs.mx <- gs.mx * scale.div
  gs.mx[which(gs.mx > clamp)] <- clamp
  gs.mx <- gs.mx / max(gs.mx)
  
  ## Plot normalized diversity index values in embedded space.
  setwd("../")
  for(distName in distNames) {
    plotData(distName, annotation = gs.mx[, distName], colorScale = "diversity", 
             plot.info = c(plot.dir, paste0(div.type, "-K", K, ifelse(to.log, "-log", ""), 
                                            ifelse(scale.div != 1, paste0("-x", scale.div), ""),
                                            ifelse(clamp != 1, paste0("-clamp-", clamp), ""))))
  }
}

################################ Inputs #################################
# 1. Define channels of interest.
# Using everything selected for viSNE.

# 2. Select FCS file with tSNE channels.
setwd(sourceDir)
fcs.name <- list.files(pattern = ".fcs")

# 3. Read in data from FCS file (apply arsinh transformation to all 
# channels except tSNE1 and tSNE2).
fcs.data <- get.fcs.data(fcs.name, 
                         arcsinh = ifelse(data.set == "68K_PBMC" ||
                                            data.set == "Synth_10_50" ||
                                            data.set == "Synth_bifurcation", FALSE, TRUE))

# 4. Import annotation.
if (data.set == "Synth_10_50" || data.set == "Synth_bifurcation" || data.set == "Panorama") {
  annotation <- fcs.data$ComponentID
  fcs.data <- fcs.data[, -which("ComponentID" == colnames(fcs.data))]
  if (data.set == "Synth_10_50") annotation <- factor(annotation)
  if (data.set == "Panorama") annotation <- factor(annotation, levels = 2:25)
} else if (data.set == "Wanderlust") {
  annotation <- fcs.data$wanderlust
  fcs.data <- fcs.data[, -which("wanderlust" == colnames(fcs.data))]
  annotation <- round(annotation * 1000)
} else if (data.set == "68K_PBMC") {
  setwd("clustered")
  fcs.name <- list.files(pattern = ".fcs")
  annotation <- get.fcs.data(fcs.name, arcsinh = FALSE, xshift = TRUE)
  annotation <- factor(annotation$clusterID)
}

############################# Process Data ##############################
# 1. Set path to the output directory.
setwd(outputDir)

# 2. Calculate distance in tSNE space.
cells.df <- as.matrix(fcs.data[, c("tSNE1", "tSNE2")])
getDistMx(cells.df, distName = "tSNE")

# 3. Calculate distance in original space.
cells.df <- as.matrix(fcs.data[, -which(colnames(fcs.data) %in% c("tSNE1", "tSNE2"))])
getDistMx(cells.df, distName = "Original")

# 4. Calculate distance in PCA space.
mypca <- princomp(cells.df)
dist.df <- as.matrix(mypca$scores)
write.csv(dist.df, "PCA-Cell-Coordinates.csv")
## Plot PCA scree.
setEPS()
postscript("PCA-Scree-Plot.eps", height = 4, width = 4)
plot(mypca)
dev.off()
## Export PCA embedding.
cells.df <- dist.df[, 1:2]
colnames(cells.df) <- c("PC1", "PC2")
getDistMx(cells.df, distName = "PCA")

# 5. Calculate distance in scFDL space from Vortex (stretch factor = 1).
setwd(sourceDir)
cells.df <- read.csv(list.files(pattern = "scFDL"), sep = ";")
cells.df <- cells.df[order(cells.df$Index_In_File), ]
cells.df <- cells.df[, c("X", "Y")]
cells.df <- (cells.df - mean(as.matrix(cells.df))) / sd(as.matrix(cells.df))
setwd(outputDir)
write.csv(cells.df, "scFDL-Cell-Coordinates.csv")
getDistMx(cells.df, distName = "scFDL")

# 6. Calculate distance in diffusion space.
cells.df <- as.matrix(fcs.data[, -which(colnames(fcs.data) %in% c("tSNE1", "tSNE2"))])
getDiffusion(cells.df)
## Diffusion components 1-2 (DC1-DC2).
cells.df <- read.csv("Diffusion-Cell-Coordinates.csv")[, -1]
cells.df <- cells.df[, c("DC1", "DC2")]
cells.df <- (cells.df - mean(as.matrix(cells.df))) / sd(as.matrix(cells.df))
write.csv(cells.df, "DC1-DC2-Cell-Coordinates.csv")
getDistMx(cells.df, distName = "DC1-DC2")

# 7. Import distance matrix in original space.
cat("Reading in distance matrix in original space ... ")
orig.df <- read.csv(paste0("Original-", kDist, "-Dist-Matrix.csv"))[, -1]
cat(paste("Done!", "\n"))

# 8. Calculate distance in classical MDS space.
# cells.df <- cmdscale(as.dist(orig.df), eig=TRUE, k=2)$points[, 1:2]  # k is the number of dim
# colnames(cells.df) <- c("MDS1", "MDS2")
# cells.df <- (cells.df - mean(as.matrix(cells.df))) / sd(as.matrix(cells.df))
# write.csv(cells.df, "MDS-Cell-Coordinates.csv")
# getDistMx(cells.df, distName = "MDS")

# 9. Calculate distance in non-classical MDS (isoMDS = isomap MASS library).
# cells.df <- isoMDS(as.dist(orig.df), k=2)$points[, 1:2]  # k is the number of dim
# colnames(cells.df) <- c("IsoMDS1", "IsoMDS2")
# cells.df <- (cells.df - mean(as.matrix(cells.df))) / sd(as.matrix(cells.df))
# write.csv(cells.df, "IsoMDS-Cell-Coordinates.csv")
# getDistMx(cells.df, distName = "IsoMDS")

# 10. Calculate distance in non-classical MDS (isomap = isomap vegan library).
# cells.df <- isomap(as.dist(orig.df), ndim=2, k=nrow(orig.df) - 1)$points[, 1:2]  # k is the number of dissimilarities retained/cell
# colnames(cells.df) <- c("Isomap1", "Isomap2")
# cells.df <- (cells.df - mean(as.matrix(cells.df))) / sd(as.matrix(cells.df))
# write.csv(cells.df, "Isomap-Cell-Coordinates.csv")
# getDistMx(cells.df, distName = "Isomap")

# 12. Plot data in embedded space.
cat(paste("Plotting data in embedded space", "\n"))
setwd(outputDir)
plotData(distName = "tSNE", annotation = annotation)
plotData(distName = "PCA", annotation = annotation)
plotData(distName = "scFDL", annotation = annotation)
plotData(distName = "DC1-DC2", annotation = annotation)

# 13. Plot correlation between original distance and embedded distance.
cat(paste("Calculating embedding quality metrics", "\n"))
Ks <- c(100, 300, 1000, 3000, 10000, 49999)
if (nrow(orig.df) <= 10000) Ks <- c(30, 100, 300, 1000, 3000, 9999)
if (data.set == "Wanderlust" && sample.set == "all_cells") {
  Ks <- c(30, 100, 300, 1000, 3000, 10863)
}
if (nrow(orig.df) <= 1000) Ks <- c(3, 10, 30, 100, 300, 999)
getCor(orig.df, distName = "PCA", K = Ks)
getCor(orig.df, distName = "tSNE", K = Ks)
getCor(orig.df, distName = "scFDL", K = Ks)
getCor(orig.df, distName = "DC1-DC2", K = Ks)

# 14. Plot correlation results.
cat(paste("Plotting embedding quality metrics", "\n"))
setwd(outputDir)
for(type in c("dist-pearson-cor", "dist-spearman-cor",
              "rank-fraction", "knn-density")) {
  ## Read in original vs. embedded comparison results.
  cor.df <- rbind(readCor(distName = "PCA", type = type),
                  readCor(distName = "tSNE", type = type),
                  readCor(distName = "scFDL", type = type),
                  readCor(distName = "Diffusion", type = type))
  
  ## Write out resupts for all methods.
  write.csv(cor.df, paste0("Embedded-vs-Original-", type, ".csv"))
  
  ## Plot results as a violin and boxplots.
  if (type != "knn-density") plotCor(cor.df, type = type, violin = FALSE)
}

# 15. Plot results over KNN density.
cat(paste("Plotting embedding quality over KNN density buckets", "\n"))
setwd(outputDir)
dir.create("KNN-Density-Plots")
for(type in c("dist-pearson-cor", "dist-spearman-cor", "rank-fraction")) {
  plotKnnDensity(type) 
}

# 16. Plot trajectory continuity or population purity.
if (!is.factor(annotation)) {
  cat(paste("Plotting trajectory continuity", "\n"))
  setwd(outputDir)
  dir.create("Trajectory-Continuity")
  invisible(lapply(c("Original", "tSNE", "PCA", "scFDL", "Diffusion"), 
         function(x) plotTrajContinuity(distName = x, annotation)))
} else {
  cat(paste("Plotting population purity", "\n"))
  setwd(outputDir)
  dir.create("Population-Purity")
  distNames <- c("tSNE", "PCA", "scFDL", "DC1-DC2")
  #divTypes <- c("Gini-Simpson", "1-Shannon", "Inv-Simpson")
  divTypes <- "Gini-Simpson"
  invisible(lapply(distNames, function(x) getPopPurity(distName = x, annotation, Ks)))
  for (div.type in divTypes) {
    invisible(lapply(Ks[1:3], function(x) plotPopPurity(distNames, annotation, K = x, 
                                                   div.type = div.type, clamp = 1)))
    invisible(lapply(Ks[1:3], function(x) plotPopPurity(distNames, annotation, K = x, 
                                                   div.type = div.type, clamp = 0.75)))
    invisible(lapply(Ks[1:3], function(x) plotPopPurity(distNames, annotation, K = x, 
                                                   div.type = div.type, clamp = 0.5)))
  }
}
