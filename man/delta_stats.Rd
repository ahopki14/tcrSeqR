\name{delta_stats}
\alias{delta_stats}
\title{Calculate change in metrics}
\description{Calculates the change in metrics between two sample types, for all patients (for 
example) in an immunoSeqR plot_ds object.}
\usage{out <- delta_stats(plot_ds,s1,s2,col,metrics,split_on)}
\arguments{
  \item{plot_ds}{An immunoSeqR plot_ds object}
  \item{s1}{The name of the factor level for the first sample to be compared (e.g. 'PRE')}
  \item{s2}{The second sample to be compared (e.g. 'POST')}
  \item{col}{Which column of plot_ds to look in for s1 and s2}
  \item{metrics}{A list of metrics (column names in plot_ds) to be compared}
  \item{split_on}{The column of plot_ds containing the objects to be looped over (e.g. 'patient')}

}
\details{}
\value{
  \item{out}{A data frame containing rows for every object in <split_on> and columns for every 
object in <metrics>. The values are log2(s2/s1).}
}
\references{}
\author{
Alexander Hopkins
}
\note{}

\seealso{}
\examples{}
\keyword{}
