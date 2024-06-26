# triplot.rda() function
#
# A function to draw a triplot (scaling 1 or scaling 2) from an object 
# of class "rda" containing RDA result from vegan's rda() function.
# 
# This version can handle RDA results with a single or more canonical axes. It can also
# draw any combination of canonical and residual axes, or pairs of canonical axes or pairs
# of residual axes.
# This function requires package vegan.
#
# ARGUMENTS
#
# ##### General parameters
# res.rda          An rda{vegan} object.
# ax1, ax2         Canonical axes to be drawn as abscissa and ordinate. Defaults: 1 and 2.
# site.sc          Can be set to "lc" (linear constraints or model scores, default) 
#                  or "wa" (weighted averages, default in vegan).
# scaling          Scaling type: only 1 or 2 are supported. Default: 1.
#
# ##### Items to be plotted
# plot.sites       If TRUE, the sites will be plotted as small circles.
# plot.spe         If TRUE, the species (or other response variables) will be plotted.
# plot.env         If TRUE, arrows for the explanatory variables will be plotted.
# plot.centr       If TRUE, symbols are plotted at the centroids of factor levels, if any. 
# arrows.only      if TRUE, plot arrows for quant. explanatory var. and factor classes
# label.sites      If TRUE, labels are added to the site symbols.
# label.spe        If TRUE, labels are added to the species arrows.
# label.env        If TRUE, labels are added to the environmental variable arrows.
# label.centr      If TRUE, labels are added to the centroids of factor levels.
# cex.char1        Character size (for sites and response variables).
# cex.char2        Character size (for explanatory variables).
# cex.point        Size of points representing centroids of factor levels.
# pca.ax           If there is a single canonical axis, use a PCA axis of the residuals
#                  as the ordinate of the triplot. Default: pca.ax=1
#
# ##### Label positions, colour centroids
# ## Positions: 1=below the point, 2=left, 3=above, 4=right. Default: 4.
# ## Note - Argument pos=NULL centres the label on the position of the object (site point,  
# ## species or environmental variable arrow, centroid) when the object is not drawn.
# pos.sites        Position of site labels. 1 to 4, as above. Default: 2.
# pos.spe          Position of species labels. 1 to 4, as above. Default: 4.
# pos.env          Position of env.variable labels. 1 to 4, as above. Default: 4.
# pos.centr        Position of centroid labels. 1 to 4, as above. Default: 4. 
# col.centr        Colour of centroid labels. Default: "blue". 
#                  Change to another colour as needed, for ex. "turquoise4"
#
# ##### Multipliers, selection of species to be plotted
# mult.spe         Multiplier for length of the species arrows. Default: 1.
# mult.arrow       Multiplier for length of the environmental arrows. Default: 1.
# select.spe       Vector containing a selection of the species numbers to be drawn in 
#                  the biplot, e.g. c(1,2,5,8,12). Draw all species if select.spe=NULL 
#                  (default value). The species that are well represented in the RDA plot 
#                  can be identified using goodness(RDA.output.object,display="species")
#
# ##### Position of the plot in frame, margins
# mar.percent      Factor to expand plot size to accomodate all items and labels. Positive 
#                  values, 0.g. 0.1, increase the margins around the plot, negative  
#                  values, e.g. -0.1, reduce them.
# optimum          If TRUE, the longest species and environmental arrows are stretched to 
#                  a length equal to the distance to the origin of the site farthest from 
#                  the origin in the plot of (ax1,ax2). This is an optimal combined 
#                  representation of the three elements. The lengths of the species and 
#                  environmental arrows can be further modified using the arguments 
#                  mult.spe and mult.arrow.
# move.origin      Move plot origin right-left and up-down. Default: move.origin=c(0,0).
#                  Ex. move.origin=c(-1,0.5) moves origin by 1 unit left and 0.5 unit up.
#
# ##### Varia
# silent           If FALSE, intermediate computation steps are printed. Default: TRUE.
#
# 
# # Example 1 - Table 11.3 of Legendre & Legendre (2012, p. 644), first 6 species only
#
# Y.mat = matrix(c(1,0,0,11,11,9,9,7,7,5,0,0,1,4,5,6,7,8,9,10,0,0,0,0,17,0,13,0,10,0,0, 
# 0,0,0,7,0,10,0,13,0,0,0,0,8,0,6,0,4,0,2,0,0,0,1,0,2,0,3,0,4),10,6)
# rownames(Y.mat) = paste("site",1:10,sep="")
# colnames(Y.mat) = paste("Sp",1:6,sep="")
# Depth = 1:10
# Sub. = as.factor(c(rep(1,3),3,2,3,2,3,2,3))     # 1=Sand, 2=Coral, 3=Other
# env = cbind(data.frame(Depth),data.frame(Sub.))
# 
# rda.out = rda(Y.mat~ Depth+Sub.,env)
# 
# # Scaling=1
# par(mfrow=c(1,2))
# triplot.rda(rda.out, scaling=1, mar.percent=0)
# triplot.rda(rda.out, scaling=1, move.origin=c(5,-5), mar.percent=-0.1)
#
# # Scaling=2
# par(mfrow=c(1,2))
# triplot.rda(rda.out, scaling=2, mar.percent=0.15, silent=FALSE)
# triplot.rda(rda.out, scaling=2, move.origin=c(0.4,-0.25), mar.percent=0.05,silent=FALSE)
# With silent=FALSE, additional information is provided on the axes positions and ranges
#
# # Example 2 - Same data
# 
# # A single quantitative explanatory variable; a single canonical axis is produced
# ( rda.out.1 = rda(Y.mat~ Depth, data=env) )
# 
# # Plot lc-scores. Compare vegan's plot.cca() and triplot.rda()
# par(mfrow = c(1, 2))
# plot(rda.out.1, scaling=1, display=c("sp","lc","cn"), col = "green4", cex = .7) 
# triplot.rda(rda.out.1, scaling=1)
# 
# # Any of the residual PCA axes can be used as the ordinate. Recommended: use ax2 = 2
#   i.e., the first non-canonical axis.
# par(mfrow = c(1, 2))
# triplot.rda(rda.out.1, scaling=1, ax1 = 1, ax2 = 2) # Ordinate: first residual axis
# triplot.rda(rda.out.1, scaling=2, ax1 = 1, ax2 = 3) # Ordinate: second residual axis
# 
# # Example 3 - Dune data from vegan
# 
# library(vegan)
# data(dune)
# data(dune.env)
# 
# rda.dune = rda(dune ~ .,dune.env)
# 
# Select the species well represented on RDA axes 1 and 2 of the triplot
# tmp = goodness(rda.dune)
# ( sp.sel = which(tmp[,2] >= 0.4) )
#
# # Scaling=2. Left: plot witl all species. Right: with selected species
# par(mfrow=c(1,2))
# triplot.rda(rda.dune, scaling=2, mar.percent=0)
# triplot.rda(rda.dune, scaling=2, 
#       select.spe=sp.sel, move.origin=c(-0.3,0), mar.percent=0.1)
#
# Use all species but plot the two first residual (non-canonical) axes, i.e., axes 13 
# and 14. Left plot with scaling = 1, right plot with scaling 2.
# dev.new(noRStudioGD = TRUE, width = 12, height = 6)
# par(mfrow = c(1,2))
# triplot.rda(rda.dune, ax1 = 13, ax2 = 14, scaling = 1)
# triplot.rda(rda.dune, ax1 = 13, ax2 = 14, scaling = 2)
#
######################################################################################
# Users interested in similar plots but in ggplot-style may try the sister function 
# ggtriplotRDA.R. This function requires the additional packages tidyverse and ggrepel
# with their dependencies.
######################################################################################
#
# #####
#
# License: GPL-2 
# Authors: Francois Gillet, Daniel Borcard & Pierre Legendre, 2016, 2021
#          This version is dev5.0.1 (4 June 2021)

'triplot.rda' <-
   function(res.rda, 
      ax1=1,
      ax2=2,
      site.sc="lc",
      scaling=1,
      plot.sites=TRUE, 
      plot.spe=TRUE,
      plot.env=TRUE,
      plot.centr=TRUE, 
      arrows.only=FALSE,
      label.sites=TRUE,
      label.spe=TRUE, 
      label.env=TRUE, 
      label.centr=TRUE,
      cex.char1=0.7, 
      cex.char2=1.0,
      cex.point=0.8,
      pos.sites=2,
      pos.spe=4,
      pos.env=4,  
      pos.centr=4,
      col.centr="blue",
      mult.spe=1, 
      mult.arrow=1,
      select.spe=NULL,
      mar.percent=0.15,  
      optimum=TRUE,
      move.origin=c(0,0),
      silent=TRUE) {

### Internal functions
#
'stretch' <- 
   function(sites, mat, ax1, ax2, n, silent=silent) {
   # Compute stretching factor for the species or environmental arrows
   # First, compute the longest distance to centroid for the sites
   tmp1 <- rbind(c(0,0), sites[,c(ax1,ax2)])
   D <- dist(tmp1)
   target <- max(D[1:n])
   # Then, compute the longest distance to centroid for the species or environmental arrows
   if (is.matrix(mat)) {
     p <- nrow(mat)   # Number of species or env. arrows to be drawn
     tmp2 <- rbind(c(0,0), mat[,c(ax1,ax2)])
     D <- dist(tmp2)
     longest <- max(D[1:p])
   } else {
     tmp2 <- rbind(c(0,0), mat[c(ax1,ax2)]) 
     longest <- dist(tmp2)
     # print(tmp2)
   }  # If a single row left in 'mat'
  #
  if(!silent)
     cat("target =", target,
         " longest =", longest,
         " fact =",target/longest,
         "\n")
  fact <- target/longest
}
#
'larger.plot' <- 
   function(sit.sc, 
            spe.sc,
            BP.sc,
            percent,
            move.origin,
            ax1,
            ax2,
            k) {
   # Internal function to expand plot limits (adapted from code by Pierre Legendre)
  if(ax1 > k) {
   mat <- rbind(sit.sc, spe.sc)
   }
  else {
   mat <- rbind(sit.sc, spe.sc, BP.sc)
   }
   range.mat <- apply(mat, 2, range)
   rownames(range.mat) <- c("Min","Max")
   z <- apply(range.mat, 2, function(x) x[2]-x[1])
   range.mat[1,] <- range.mat[1,]-z*percent
   range.mat[2,] <- range.mat[2,]+z*percent
   if(move.origin[1] != 0) range.mat[,ax1] <- range.mat[,ax1] - move.origin[1]
   if(move.origin[2] != 0) range.mat[,ax2] <- range.mat[,ax2] - move.origin[2]
   range.mat
  }
### End of internal functions

  if(class(res.rda)[1]!="rda" & class(res.rda)[2]!="rda")
     stop("The input file is not a vegan rda output object")
  if(length(res.rda$colsum)==1)
     stop("Function triplot.rda is not compatible with results that contain no species scores")
  if(scaling!=1 & scaling!=2) 
     stop("Function only available for scaling = 1 or 2")
  if(ax1 > ax2)
     stop("Axes must be drawn in increasing order")

  k <- length(res.rda$CCA$eig)        # number of RDA eigenvalues
  if(ax1 <= k & site.sc=="lc") {
     cat(  "-----------------------------------------------------------------------")
     cat("\nSite constraints (lc) selected. To obtain site scores that are weighted") 
     cat("\nsums of species scores (default in vegan), argument site.sc must be set")
     cat("\nto wa.")
     cat("\n-----------------------------------------------------------------------\n\n")
                    }
  n.sp <- length(res.rda$colsum)      # number of species
  # 'vec' will contain the selection of species to be drawn
  if (is.null(select.spe)) {
    vec <- 1:n.sp } else {
    vec <- select.spe
  }

 ahead <- 0.05   		# Length of arrow heads
 aangle <- 30    		# Angle of arrow heads
	
# Properties of species and site scores in scalings 1 and 2 –
# Scaling 1: the species scores have norms of 1
# Scaling 1: the site scores are scaled to variances = can.eigenvalues
# Scaling 2: the species scores have norms of sqrt(can.eigenvalues)
# Scaling 2: the site scores are scaled to variances of 1
# --------------------------------------------------------------------

### This version reconstructs the original RDA output in L&L 2012, Section 11.1.3

  Tot.var <- res.rda$tot.chi            # Total variance in response data Y
  eig.val <- c(res.rda$CCA$eig, res.rda$CA$eig) # all eigenvalues
  Lambda <- diag(eig.val)               # Diagonal matrix of eigenvalues
  eig.val.rel <- eig.val / Tot.var      # Relative eigenvalues of Y-hat
  Diag <- diag(sqrt(eig.val.rel))       # Diagonal matrix of sqrt(relative eigenvalues)
  U.sc1 <- cbind(res.rda$CCA$v, res.rda$CA$v) # All species scores, scaling=1
  U.sc2 <- U.sc1 %*% sqrt(Lambda)       # Species scores, scaling=2
  colnames(U.sc2) <- colnames(U.sc1)
  n <- nrow(res.rda$CCA$u)              # Number of observations
  Z.sc2 <- cbind(res.rda$CCA$u, res.rda$CA$u) * sqrt(n - 1)  # "lc" site scores, scaling=2
  Z.sc1 <- Z.sc2 %*% sqrt(Lambda)       # "lc" site scores, scaling=1
  colnames(Z.sc1) <- colnames(Z.sc2)
  F.sc2 <- cbind(res.rda$CCA$wa, res.rda$CA$u) * sqrt(n - 1)  # "wa" site scores, scaling=2
  F.sc1 <- F.sc2 %*% sqrt(Lambda)       # "wa" site scores, scaling=1
  colnames(F.sc1) <- colnames(F.sc2)
  BP.sc2 <- res.rda$CCA$biplot          # Biplot scores, scaling=2 ; cor(Z.sc1, X)
  BP.sc2 <- cbind(BP.sc2, matrix(0, nrow = nrow(BP.sc2), 
                                 ncol = length(eig.val) - k))
  colnames(BP.sc2) <- colnames(F.sc2)
  BP.sc1 <- BP.sc2 %*% Diag             # Biplot scores, scaling=1
  colnames(BP.sc1) <- colnames(BP.sc2)

  if (!is.null(res.rda$CCA$centroids)) {
    centroids.sc2 <- res.rda$CCA$centroids * sqrt(n - 1) # Centroids, scaling=2
    centroids.sc2 <- cbind(centroids.sc2, matrix(0, nrow = nrow(centroids.sc2), 
                                                 ncol = length(eig.val) - k))
    colnames(centroids.sc2) <- colnames(F.sc2)
    centroids.sc1 <- centroids.sc2 %*% sqrt(Lambda)      # Centroids, scaling=1
    colnames(centroids.sc1) <- colnames(centroids.sc2)
  }
  centroids.present <- TRUE
  if (is.null(res.rda$CCA$centroids)) {
    centroids.present <- FALSE
    if (plot.centr | label.centr) {
      cat("\nNo factor, hence levels cannot be plotted with symbols;")
      cat("\n'plot.centr' is set to FALSE\n")
      plot.centr  <- FALSE
      label.centr <- FALSE
    }
  }
  
  if (is.null(select.spe)) {vec <- 1:n.sp} else {vec <- select.spe}
  
  if (scaling == 1) {
    if (site.sc == "lc") {
      sit.sc <- Z.sc1
    } else {
      sit.sc <- F.sc1
    }
    spe.sc <- U.sc1[vec, ]
    BP.sc  <- BP.sc1
    if (centroids.present)
      centroids <- centroids.sc1
  } else {
    # For scaling 2
    if (site.sc == "lc") {
      sit.sc <- Z.sc2
    } else {
      sit.sc <- F.sc2
    }
    spe.sc <- U.sc2[vec, ]
    BP.sc  <- BP.sc2
    if (centroids.present)
      centroids <- centroids.sc2
  }
  
  fact.spe <- 1
  fact.env <- 1
  if (centroids.present & (plot.centr | label.centr)) {
    to.plot <- which(!(rownames(BP.sc) %in% rownames(centroids)))
  } else {
    to.plot <- 1:nrow(BP.sc)
  }
  
  if (optimum) {
    fact.spe <-
      stretch(sit.sc, spe.sc, ax1, ax2, n, silent = silent)
    if (arrows.only) {
      fact.env <-
        stretch(sit.sc, BP.sc, ax1, ax2, n, silent = silent)
    } else {
      # arrows only==FALSE
      quant.env.present <- FALSE
      if (length(to.plot) > 0) {
        quant.env.present <- TRUE
        fact.env <-
          stretch(sit.sc, BP.sc[to.plot, ], ax1, ax2, n, silent = silent)
      }
    }
  }
  
  if (!silent)
    cat("fac.spe =", fact.spe, "   fact.env =", fact.env, "\n")
  spe.sc <- spe.sc * fact.spe * mult.spe
  BP.sc <- BP.sc * fact.env * mult.arrow
  lim <-
    larger.plot(
      sit.sc[, ],
      spe.sc[, ],
      BP.sc[, ],
      percent = mar.percent,
      move.origin = move.origin,
      ax1 = ax1,
      ax2 = ax2,
      k = k
    )
  if (!silent)
    print(lim)


### Drawing the triplot begins here ###

# Draw the main plot
  mat <- rbind.data.frame(sit.sc, spe.sc, BP.sc)
#  mat <- rbind(sit.sc[,1:k], spe.sc[,1:k], BP.sc[,1:k])

  if(ax1 <= k){
    titre <- "RDA triplot - Scaling"
              }
  else{
    titre <- "Biplot of residuals of RDA - Scaling"
      }    

  plot(mat[,c(ax1,ax2)], 
    type="n", 
    main=paste(titre, scaling, "-", site.sc),
    xlim=c(lim[1,ax1], lim[2,ax1]),
    ylim=c(lim[1,ax2], lim[2,ax2]), 
    xlab=colnames(U.sc1)[ax1], 
    ylab=colnames(U.sc1)[ax2],
    asp=1)
  abline(h=0, 
    v=0, 
    col="grey60")
  
  # Draw the site scores ("lc" or "wa" for RDA, implicitly "wa" only for residuals (PCA))
  if(plot.sites) {
    points(sit.sc[,ax1], sit.sc[,ax2], pch=20)
    if(label.sites)
      text(sit.sc[,ax1], sit.sc[,ax2], labels = rownames(sit.sc), col="black", pos=pos.sites, cex=cex.char1)
  } else {
    if(label.sites)
      text(sit.sc[,ax1], sit.sc[,ax2], labels = rownames(sit.sc), col="black", pos=NULL, cex=cex.char1)
  }

  # Draw the species scores
  if(plot.spe) {
    arrows(0, 0, spe.sc[,ax1], spe.sc[,ax2], length=ahead, angle=aangle, col="red")
  	if(label.spe)
    text(spe.sc[,ax1], spe.sc[,ax2], labels = rownames(spe.sc), col="red", pos=pos.spe, cex=cex.char1)
  } else {
  	if(label.spe)
    text(spe.sc[,ax1], spe.sc[,ax2], labels = rownames(spe.sc), col="red", pos=NULL, cex=cex.char1)
  }

  # Draw the explanatory variables (only if at least one axis is canonical)
  #
  if(ax1 <= k) {
  if(!arrows.only) {
  # 1. Quantitative variables
    if(quant.env.present & plot.env) {   # Print arrows and labels for quantitative var.
      arrows(0, 0, BP.sc[to.plot,ax1]*mult.arrow, BP.sc[to.plot,ax2]*mult.arrow, length=ahead, angle=aangle, col="blue")
  	  if(label.env)   # Print labels for the quantitative variables
        text(BP.sc[to.plot,ax1]*mult.arrow, BP.sc[to.plot,ax2]*mult.arrow, labels = rownames(BP.sc)[to.plot], col="blue", pos=pos.env, cex=cex.char2)
    } else {
  	if(quant.env.present & !plot.env & label.env)   # Only print labels for quant. var.
        text(BP.sc[to.plot,ax1]*mult.arrow, BP.sc[to.plot,ax2]*mult.arrow, labels = rownames(BP.sc)[to.plot], col="blue", pos=NULL, cex=cex.char2)
    }
  #
  # 2. Centroids and labels of factor levels
    if(centroids.present & plot.centr) {   # Print symbols and labels for factor classes
      points(centroids[,ax1], centroids[,ax2], pch=18, cex=cex.point, col=col.centr)
      if(label.centr)
        text(centroids[,ax1], centroids[,ax2], labels = rownames(centroids), col="blue", pos=pos.centr, cex=cex.char2)
    } else {
    #
      if(centroids.present & !plot.centr & label.centr)   # Only print labels for classes
        text(centroids[,ax1], centroids[,ax2], labels = rownames(centroids), col="blue", pos=NULL, cex=cex.char2)
    }
  }

  # 3. All env. var.: plot arrows and labels for all var. in 'BP.sc', quant. and factors
  if(arrows.only) {
    arrows(0, 0, BP.sc[,ax1]*mult.arrow, BP.sc[,ax2]*mult.arrow, length=ahead, angle=aangle, col="blue")
  	if(label.env)  # Print labels for the quantitative variables
      text(BP.sc[,ax1]*mult.arrow, BP.sc[,ax2]*mult.arrow, labels = rownames(BP.sc), col="blue", pos=pos.env, cex=cex.char2)
  }
  }
#
}
