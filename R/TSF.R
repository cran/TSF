#=====================================================================================#
# PURPOSE : Application 0f Two Stage Forecasting Approach in Long Memory Time Series  #
# AUTHOR  : Sandipan Samanta, Ranjit Kumar Paul and Dipankar Mitra                    #
# DATE    : 14 July, 2017                                                             #
# VERSION : Ver 0.1.0                                                                 #
#=====================================================================================#

fdseries <- function(x, d)
{
  x <- as.data.frame(x)
  names(x) <- "series"
  x <- x$series
  if (NCOL(x) > 1)
    stop("only implemented for univariate time series")
  if (any(is.na(x)))
    stop("NAs in x")
  n <- length(x)
  stopifnot(n >= 2)
  PI <- numeric(n)
  PI[1] <- d
  for (k in 2:n) {
    PI[k] <- PI[k-1]*(d - k + 1)/k
  }
  FractionalDiffSeries <- x
  for (i in 2:n) {
    FractionalDiffSeries[i] <- x[i] + sum(PI[1:(i-1)]*x[(i-1):1])
  }
  return(FractionalDiffSeries)
}
StructuralBrekwithLongmemory <- function(ts,bandwidth)
{
  d_ <- fdGPH(ts, bandw.exp = bandwidth)$d #fracdiff package is used for GPH estimate of long memory

  r=fdseries(ts, d=d_) # to obtain fractional differenced series, diffseries function defined above is used
  f=auto.arima(r)# to obtain one step ahead forecast of fractional differenced series, auto.arima is used

  fpred <- fitted(f)

  ff=forecast(f,h=1)
  #' Returing to orginal series
  updatedSeries <- fdseries(rbind(as.matrix(r),as.matrix(ff$mean[1])),d=-d_)
  predictedSeries <- fdseries(rbind(as.matrix(fpred)),d=-d_)
  return(list(updatedSeries=updatedSeries,predictedSeries=predictedSeries,LongMemoryParam=d_))
}
forecastTSF <- function(N0,Xt,bandwidth)
{
  OutSBwLM <- StructuralBrekwithLongmemory(Xt,bandwidth)
  updatedts <- OutSBwLM$updatedSeries; Prediction <- OutSBwLM$predictedSeries; AllLongMemoryParam <- OutSBwLM$LongMemoryParam;
  for(i in 2 : N0){
    OutSBwLM <- StructuralBrekwithLongmemory(updatedts,bandwidth)
    updatedts <- OutSBwLM$updatedSeries; updatedLongmemory <- OutSBwLM$LongMemoryParam;
    forecastseries<-updatedts[(length(updatedts) - N0 + 1):length(updatedts)]
    AllLongMemoryParam <- c(AllLongMemoryParam,updatedLongmemory)
  }
  return(list(forecastseries=forecastseries,Prediction=Prediction,LongMemoryParameter=AllLongMemoryParam))

}
