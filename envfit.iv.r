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
