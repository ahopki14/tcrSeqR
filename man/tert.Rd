\name{tert}
\alias{tert}
\title{Make tertiles}
\description{A function to make tertiles from a vector}
\usage{tert <- tert(x)}
\arguments{
  \item{x}{A numeric vector}
}
\details{}
\value{
  \item{tert}{A vector of the same length as x, containing factors 1, 2 and 3, indicating in 
which tertile the data belongs.}
}
\references{}
\author{
Alexander Hopkins
}
\note{}

\seealso{}
\examples{
set.seed(1)
x <- runif(12)
plot(tert(x),x)
}
\keyword{}
