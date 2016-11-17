# Requires quantreg.
library(quantreg)

# Robust fitting with the in-built R line function.
tukeyline <- function(y,k=5) {
  m <- length(y)
  x <- -k:k
  r <- vector("numeric",m)
  for (i in (k+1):(m-k)) {
    z <- y[(i - k):(i + k)]
    r[i] <- ifelse(sum(!is.na(z))>2,line(x=x,y=z)$fitted.values[k+1],NA)
  }
  return(r)
}

# Quantile (median) regression with the rq function of the quantreg R package.
medianreg <- function(y,k=5) {
  m <- length(y)
  x <- -k:k
  r <- vector("numeric",m)
  for (i in (k+1):(m-k)) {
    z <- y[(i - k):(i + k)]
    r[i] <- ifelse(sum(!is.na(z))>5,rq(z~x,.5)$fitted.values[k+1],NA)
  }
  return(r)
}

runmedian <- function(y,k=5) {
  n <- length(y)
  return(c(y[1:k],sapply((k+1):(n-k),function(j) median(y[(j - k):(j + k)],na.rm=T)),y[(n-k+1):n]))
}

# Hampel filter, modified to not break down when NA are present in the data.
# Options include the possibility of using the tukeyline or the quantile regression
# estimations, instead of the local median.

hampel.filter <- function(y,k=5,t0=3,method="hampel",reverse=T) {
  n <- length(y)
  if (reverse) {
    z <- c(rev(y[2:(k+1)]),y,rev(y[(n-k):(n-1)])) 
    n <- length(z)
  } else z <- y
  z0 <- switch(method,
               hampel=runmedian(y=z,k=k),
               tukey=tukeyline(y=z,k=k),
               quantile=medianreg(y=z,k=k))
  z0.na <- !is.na(z0)
  z.na <- !is.na(z)
  ind <- NULL
  L <- 1.4826
  for (i in (k + 1):(n - k)) {
    S0 <- L * median(abs(z[(i - k):(i + k)] - z0[i]),na.rm=T)
    if (z0.na[i] & z.na[i] & !is.na(S0) & abs(z[i] - z0[i]) > t0 * S0) {
      z[i] <- z0[i]
      ind <- c(ind, i)
    }
  }
  if (reverse) {
    z <- z[(k+1):(n-k)]
    ind <- ind-k
  }
  return(list(y = z, ind = ind))
}