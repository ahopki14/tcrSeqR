\name{mean_95}
\alias{mean_95}
\title{Computes mean and 95\% CI}
\description{To be used with stat_summary to make error bars centered at the mean, and 
covering the 95\% confidence interval.}
\usage{m <- mean_95(x)}
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

\seealso{mean_sd}
\examples{
ggplot(...) + geom_point(...) + stat_summary(fun.data=mean_95)
}
\keyword{}
