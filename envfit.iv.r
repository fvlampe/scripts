envfit.iv <- function (ord, veg, spec.iv, permutations = 999, choices = c(1, 2), display = "sites", w = weights(ord), na.rm = FALSE, ...) 
{
  weights.default <- function(object, ...) NULL
  vectorfit.iv <-
    function (X, veg, spec.iv, permutations, w, ...)
    {
      apply.FUN <- function (x) 
      {
        veg.temp <- veg [,!is.na (x)]
        x.temp <- x[!is.na (x)]
        colSums (t(veg.temp)*x.temp)/rowSums (veg.temp)
      }
      
      apply.FUN.sample <- function (x) 
      {
        veg.temp <- veg [,!is.na (x)]
        x.temp <- x[!is.na (x)]
        colSums (t(veg.temp)*sample (x.temp))/rowSums (veg.temp)
      }
      
      P <- apply (spec.iv, 2, FUN = apply.FUN)
      X <- as.matrix(X)
      if (missing(w) || is.null(w)) 
        w <- 1
      if (length(w) == 1) 
        w <- rep(w, nrow(X))
      Xw <- .Call("do_wcentre", X, w)
      dim(Xw) <- dim(X)
      Pw <- .Call("do_wcentre", P, w)
      dim(Pw) <- dim(P)
      colnames(Pw) <- colnames(P)
      nc <- ncol(X)
      Q <- qr(Xw)
      H <- qr.fitted(Q, Pw)
      heads <- qr.coef(Q, Pw)
      r <- diag(cor(H, Pw)^2)
      heads <- decostand(heads, "norm", 2)
      heads <- t(heads)
      if (is.null(colnames(X))) 
        colnames(heads) <- paste("Dim", 1:nc, sep = "")
      else colnames(heads) <- colnames(X)
      if (permutations) {
        nr <- nrow(X)
        permstore <- matrix(nrow = permutations, ncol = ncol(P))
        for (i in 1:permutations) {
          take <- apply (spec.iv, 2, FUN = apply.FUN.sample)
          take <- .Call("do_wcentre", take, w)
          dim(take) <- dim(P)
          Hperm <- qr.fitted(Q, take)
          permstore[i, ] <- diag(cor(Hperm, take))^2
        }
        permstore <- sweep(permstore, 2, r, ">")
        pvals <- (apply(permstore, 2, sum) + 1)/(permutations + 
                                                   1)
      }
      else pvals <- NULL
      sol <- list(arrows = heads, r = r, permutations = permutations, 
                  pvals = pvals)
      class(sol) <- "vectorfit"
      sol
    }
  w <- eval(w)
  vectors <- NULL
  factors <- NULL
  seed <- NULL
  X <- scores(ord, display = display, choices = choices, ...)
  keep <- complete.cases(X)
  if (any(!keep)) {
    if (!na.rm) 
      stop("missing values in data: consider na.rm = TRUE")
    X <- X[keep, , drop = FALSE]
    na.action <- structure(seq_along(keep)[!keep], class = "omit")
  }
  vectors <- vectorfit.iv(X, veg, spec.iv, permutations, choices,
                          w = w, ...)
  sol <- list(vectors = vectors, factors = factors)
  if (!is.null(na.action)) 
    sol$na.action <- na.action
  class(sol) <- "envfit"
  sol
}





ordiselect_factor <- function (matrix, ord, ablim = 0.9, fitlim = 1, choices = c(1,
                  2), method = "axes", env, p.max = 0.05, freq = FALSE,
                               na.rm = FALSE)
{
  if (!is.data.frame(matrix)) {
    matrix <- data.frame(matrix)
  }
  if (freq == F) {
    abund <- apply(matrix, 2, sum)
  }
  else {
    abund <- apply(matrix > 0, 2, sum)
  }
  if (method == "axes") {
    scores <- data.frame(scores(ord, display = "species",
                                choices = choices))
    scores1 <- scores[, 1]
    scores2 <- scores[, 2]
    rownames(scores[abund >= quantile(abund, 1 - ablim) &
                      (scores1 >= quantile(scores1, 1 - 0.5 * fitlim, na.rm = na.rm) |
                         scores1 <= quantile(scores1, 0.5 * fitlim, na.rm = na.rm) |
                         scores2 >= quantile(scores2, 1 - 0.5 * fitlim,
                                             na.rm = na.rm) | scores2 <= quantile(scores2,
                                                                                  0.5 * fitlim, na.rm = na.rm)), ])
  }
  else if (method == "vars") {
    if (class(env) != "envfit") {
      print("FATAL: Fitted environmental variables are no result of envfit()")
    }
    else {
      scores_spec <- data.frame(scores(ord, display = "species",
                                       choices = choices))
      sig <- env$factors$pvals
      scores_env <- data.frame(scores(env, display = "factors"))
      scores_env <- scores_env[sig < p.max, ]
      if (nrow(scores_env) == 0) {
        print("WARNING: No significant environmental variables. Only abundance limit used for selection.")
        names(abund[abund >= quantile(abund, 1 - ablim)])
      }
      else if (nrow(scores_env) == 1) {
        euclid <- data.frame(fields::rdist(scores_spec,
                                           scores_env))
        names(euclid) <- rownames(scores_env)
        rownames(euclid) <- names(matrix)
        rownames(data.frame(euclid[abund > quantile(abund,
                                                    1 - ablim) & (euclid <= quantile(unlist(euclid),
                                                                                     0.5 * fitlim, na.rm = na.rm) | euclid >= quantile(unlist(euclid),
                                                                                                                                       1 - 0.5 * fitlim, na.rm = na.rm)), , drop = FALSE]))
      }
      else {
        euclid <- data.frame(fields::rdist(scores_spec,
                                           scores_env))
        names(euclid) <- rownames(scores_env)
        rownames(euclid) <- names(matrix)
        best_fit <- euclid[abund > quantile(abund, 1 -
                                              ablim) & (euclid <= quantile(unlist(euclid),
                                                                           0.5 * fitlim, na.rm = na.rm) | euclid >= quantile(unlist(euclid),
                                                                                                                             1 - 0.5 * fitlim, na.rm = na.rm)), ]
        rownames(best_fit[complete.cases(best_fit), ])
      }
    }
  }
  else {
    print("FATAL: Selected Method unknown")
  }
}
