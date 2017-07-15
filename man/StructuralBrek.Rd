\name{StructuralBrekwithLongmemory}
  \alias{StructuralBrekwithLongmemory}
  \title{Predicting fractionally differenced series in presence of structural break}
\description{The function is used for prediction of long memory time series in presence of structural break
}

\usage{
StructuralBrekwithLongmemory(ts,bandwidth)}
\arguments{
  \item{ts}{univariate time series}
  \item{bandwidth}{the bandwidth used in the regression equation}
   }
\value{
    \item{StructuralBrekwithLongmemory}{the updated series at first step of TSF appraoch, prediction based on TSF approach and the estimate of long memory parameter}
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
bandwidth<-0.9
set.seed(Seed)
Sim.Series <- fracdiff::fracdiff.sim(n = N, ar = c(PHI), ma = c(THETA),
d = D, rand.gen = rnorm, sd = SD, mu = M)

Xt <- as.ts(Sim.Series$series)

## Forecasting using TSF method
StructuralBrekwithLongmemory(Xt,bandwidth)
}
