\name{iseqr_check}
\alias{iseqr_check}
\title{Check immunoSeqR ds and stats objects for alignment with the dictionary}
\description{This function checks that the dataset and, optionally the stats object are in 
the same order as the dictionary.}
\usage{check <- iseqr_check(dict,ds,stats)}
\arguments{
  \item{dict}{An immunoSeqR dictionary}
  \item{ds}{An immunoSeqR dataset}
  \item{stats}{An immunoSeqR stats object}
}
\details{}
\value{
  \item{check}{A list of two logicals, indicating if the ds and stats are in correct order.}
}
\references{}
\author{
Alexander Hopkins
}
\note{If the stat object is omitted, then check$stats will return NA}

\seealso{}
\examples{}
\keyword{}
