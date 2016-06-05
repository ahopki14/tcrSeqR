\name{iseqr_summarize}
\alias{iseqr_summarize}
\title{Performs all pairwise t-test between groups}
\description{Performs pairwise t-tests of all the groups of any metadata field, using any 
data column of an immunoSeqR plot_ds object.}
\usage{out <- iseqr_summarize(plot_ds,slpit_on,metric)}
\arguments{
  \item{plot_ds}{An immunoSeqR plot_ds object, or any data frame with metadata and data}
  \item{split_on}{The name of the column of plot_ds containing the groups to be compared}
  \item{metric}{The name of the column of plot_ds containting the data}
}
\details{}
\value{
  \item{out}{An NxN matrix, where N is the number of different elements in <split_on>, 
containing the unadjusted p-values of the comparisons.}
}
\references{}
\author{
Alexander Hopkins
}
\note{}

\seealso{}
\examples{}
\keyword{}
