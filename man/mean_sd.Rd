\name{mean_sd}
\alias{mean_sd}
\title{Computes mean and SD}
\description{To be used with stat_summary to make error bars centered at the mean, and 
covering the +/- one standard deviation.}
\usage{m <- mean_sd(x)}
\arguments{
  \item{x}{numeric vector}
}
\details{}
\value{
  \item{m}{A named numeric containing the mean, upper and lower bounds.}
}
\references{}
\author{
Alexander Hopkins
}
\note{}

\seealso{mean_95}
\examples{
ggplot(...) + geom_point(...) + stat_summary(fun.data=mean_sd)
}
\keyword{}
