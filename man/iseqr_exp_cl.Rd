\name{iseqr_exp_cl}
\alias{iseqr_exp_cl}
\title{Find Expanded Clones}
\description{Convenience function for exp_clones. Finds clones which significantly change 
between two types of samples, across multiple patients, in an immunoSeqR dataset}
\usage{out <- iseqr_exp_cl(ds,dict,s1,s2,category,by)}
\arguments{
  \item{ds}{An immunoSeqR dataset}
  \item{dict}{An immunoSeqR dictionary}
  \item{s1}{The first level of <category> to be compared} 
  \item{s2}{The second level of <category> to be compared} 
  \item{category}{The metadata field in which to look for s1 and s2} 
  \item{by}{the experimental unit to be iterated over} 

}
\details{}
\value{
  \item{out}{A list with various useful items}
  \item{out$num_exp}{A list of all the units in <by> and the number of changed clones}
  \item{out$p}{The vectors of p-values for all the comparisons}
  \item{out$big_out}{A list containing, for each unit in <by>, the counts and p-values}
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
use their code. It also returns different results from their code, so use with caution. 

Additionally, this output object can be very large, as it contains the entire dataset with 
additional p-values, etc. }

\seealso{exp_clones}
\examples{}
\keyword{}
