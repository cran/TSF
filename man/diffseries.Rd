\name{TSF}
\alias{fdseries}
\title{Fractionally differenced series for any value of d}
\description{
The function fdseries computes the fractional differenced series for any value of d i.e. positive or negetive.
}

\usage{
fdseries(x, d)}
\arguments{
  \item{x}{univariate time series}
  \item{d}{The orer of fractional differencing to be done}
   }
\value{
    \item{fdseries}{fractionally differenced series for both positive as well as negetive d}
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

set.seed(Seed)
Sim.Series <- fracdiff::fracdiff.sim(n = N, ar = c(PHI), ma = c(THETA),
d = D, rand.gen = rnorm, sd = SD, mu = M)

Xt <- as.ts(Sim.Series$series)

## fractional differencing
fdseries(Xt,d=D)
}
