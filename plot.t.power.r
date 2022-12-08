## Author: Friedemann von Lampe, 08 Dec. 2022
## Function for noncentrality parameter delta in two sample t-test 
## according to Faul et al. 2009
## d = effect size, n = sample size in groups 1 & 2
ncp.t <- function(d, n1, n2) {
  d * sqrt((n1*n2)/(n1+n2))
}



## Function for probability density distribution of H0 and H1
## in dependence of noncentrality parameter delta
## with visualisation of type I & II errors
## d = effect size, n = sample size in groups 1 & 2
plot.t.power <- function(d, n1, n2) {
  
  colalpha <- rgb(0, 0, 0, max = 255, alpha = 125, names = "gray50")
  colbeta <- rgb(255, 0, 0, max = 255, alpha = 125, names = "red50")
  
  ncp <- ncp.t(d, n1, n2)
  N <- n1 + n2
  
  if(d >= 0) {
    lo <- -6 
    up <- 6 + ncp*2         # 1/ncp*15
  } else {
    lo <- -6 + ncp*2
    up <- 6
  }
  
  xt0 <-seq(lo, up, 0.01)
  ydt0 <-dt(xt0, df = N-2)
  plot(xt0, ydt0, type = "l",
       xlab = "Value",
       ylab = "Probability",
       main = paste("t-Test power with N = ", N, "and d =", d),
       las = 1)
  
  xt1 <-seq(lo, up, 0.01)
  ydt1 <-dt(xt1, df = N-2, ncp = ncp)
  lines(xt1, ydt1, type = "l", col = "red")
  
  # Area above critical value H0
  locritvalue <-  qt(0.025, df = N-2)
  upcritvalue <-  qt(0.975, df = N-2)
  
  # Upper alpha
  polygon(c(xt0[xt0 >= upcritvalue ], upcritvalue),
          c(ydt0[xt0 >= upcritvalue ], 0),
          col = colalpha,
          border = NA)
  
  # Lower alpha
  polygon(c(xt0[xt0 <= locritvalue ], locritvalue),
          c(ydt0[xt0 <= locritvalue ], 0),
          col = colalpha,
          border = NA)
  
  # Beta
  if(d >= 0) {
    polygon(c(xt1[xt1 <= upcritvalue ], upcritvalue),
            c(ydt1[xt1 <= upcritvalue ], 0),
            col = colbeta,
            border = NA)
  }
  else {
    polygon(c(xt1[xt1 >= locritvalue ], locritvalue),
            c(ydt1[xt1 >= locritvalue ], 0),
            col = colbeta,
            border = NA)
  }
  
  abline(v = c(qt(c(0.025, 0.975), df = N-2)),
         lty = 2)
  
  legend("topright", 
         c(expression(paste(alpha," (Type I error)")), 
           expression(paste(beta," (Type II error)"))), 
         fill = c(colalpha, colbeta),
         bty = "n")
  
  legend("topleft", 
         c("H0", "Critical values", "H1"), 
         lty = c(1, 2, 1),
         col = c("black", "black", "red"),
         bty = "n")
}