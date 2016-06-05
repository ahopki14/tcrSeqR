\name{olScan}
\alias{olScan}
\title{Scans the overlap between samples}
\description{This function examines the overlap between two samples, across a variety of N, 
where N is the number of samples included when ranked by most abundant (in sample 1). 
Essentially, it first asks 'what is the overlap of the top <byN> clones?', Then the top 2 x 
<byN> clones? This is repeated until all clones are included. The result is a curve which 
approaches the actual overlap between the samples. }
\usage{out <- olScan(sample1,sample2,byN)}
\arguments{
  \item{sample1}{The first sample from an immunoSeq dataset (this is the one used for ranking) } 
  \item{sample2}{The second sample}
  \item{byN}{The step size}
}
\details{}
\value{
  \item{out}{A numeric vector, contating the overlaps}
}
\references{}
\author{
Alexander Hopkins
}
\note{
Because the first sample is used for sorting, the order of sample1 and sample2 matters.
}

\seealso{}
\examples{}
\keyword{}
