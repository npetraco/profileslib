library(fftw)

ccf.rfftw <- function(x, y, max.lag) {


  N <- max(length(x), length(y))

  #Compute the FFT size as the "next power of 2" of the input vector's length (max)
  b <- ceiling(log2(2.0 * N - 1))
  fftsize <- 2^b

  #Mean center and Pad
  xx <- c(x-mean(x), rep(0, (fftsize)-length(x) ))
  yy <- rev(c(y-mean(y), rep(0, (fftsize)-length(y) )))

  p <- planFFT(fftsize)

  #Compute CCF:
  sf <- sqrt(sum((x-mean(x)) * (x-mean(x)))) * sqrt(sum((y-mean(y)) * (y-mean(y))))
  res <- (1/sf) * Re(IFFT( FFT(xx, plan=p) * FFT(yy, plan=p), plan=p, scale=T))

  #Reformat:
  center.idx <- fftsize/2
  res <- c(res[ (length(res)/2 +1 ):length(res)], res[1:(length(res)/2)])
  res <- res[c((center.idx-max.lag):center.idx, (center.idx+1):(center.idx+max.lag))]

  lags <- seq(-max.lag,max.lag,1)

  #plot(lags, res,typ="h")

  return(list(res,lags))

}
