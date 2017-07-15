\name{forecastTSF}
  \alias{forecastTSF}
  \title{Forecasting fractionally differenced series using TSF approach}
\description{The function is used for forecasting long memory time series using TSF approach}

\usage{
forecastTSF(N0,Xt,bandwidth)}
\arguments{
  \item{N0}{lead period of forecast}
  \item{Xt}{univariate time series}
  \item{bandwidth}{the bandwidth used in the regression equation}
   }
\value{
    \item{forecastTSF}{the predicted values, the out of sample forecasts and the values of long memory parameter}
  }
\references{

Papailias, F. and Dias, G. F. 2015. Forecasting long memory series subject to structural 	change: A two-stage approach. International Journal of Forecasting, 31, 1056 to 1066.

Wang, C. S. H., Bauwens, L. and Hsiao, C. 2013. Forecasting a long memory process subject 	to structural breaks. Journal of Econometrics, 177, 171-184.

Reisen, V. A. (1994) Estimation of the fractional difference parameter in the ARFIMA(p,d,q) model using the smoothed periodogram. Journal Time Series Analysis, 15(1), 335 to 350.
}
\author{Sandipan Samanta, Ranjit Kumar Paul and Dipankar Mitra}
\examples{
## Simulating Long Memory Series
N <- 1000
PHI <- 0.2
THETA <- 0.1
SD <- 1
M <- 0
D <- 0.2
Seed <- 123
N0<-9
bandwidth<-0.9
set.seed(Seed)
Sim.Series <- fracdiff::fracdiff.sim(n = N, ar = c(PHI), ma = c(THETA),
d = D, rand.gen = rnorm, sd = SD, mu = M)

Xt <- as.ts(Sim.Series$series)

## Forecasting using TSF method
forecastTSF (N0,Xt,bandwidth)
}
