\name{olScan}
\alias{olScan}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the overlap as more samples are added
}
\description{
Calculate the overlap between two samples using progressively more samples (in descending order of abundance in the first sample)
}
\usage{
olScan(x,y,byN=10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{First data set}
  \item{y}{Second data set}
  \item{byN}{Step size}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{ol}{Vector of overlap values calculated}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Alexander Hopkins
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{

}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
