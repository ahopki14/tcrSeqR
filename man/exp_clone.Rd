\name{exp_clone}
\alias{exp_clone}
\title{Find Expanded Clones}
\description{Finds clones which significantly change between two samples in an immunoSeqR 
dataset}
\usage{p <- exp_clone(x,y)}
\arguments{
  \item{x}{The first sample (a column of ds)}
  \item{y}{The second sample } 
}
\details{}
\value{
  \item{p}{A vector of unadjusted p-values from a Fisher's Exact test, each corresponding to 
one row of x and y.}
}
\references{
DeWitt, W.S., Emerson, R.O., Lindau, P., Vignali, M., Snyder, T.M., Desmarais, C., Sanders, 
C., Utsugi, H., Warren, E.H., McElrath, J., et al. (2015). Dynamics of the Cytotoxic T Cell 
Response to a Model of Acute Viral Infection. Journal of Virology 89, 4517–4526.
}
\author{
Alexander Hopkins
}
\note{This is an implementation of the technique published by DeWitt et al., but it does not 
use their code. It also returns different results from their code, so use with caution. }

\seealso{iseqr_exp_cl}
\examples{}
\keyword{}
