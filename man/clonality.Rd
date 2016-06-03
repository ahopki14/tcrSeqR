\name{clonality}
\alias{clonality}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the clonality or entropy of a distribution
}
\description{
Calculates the clonality (1 - normalized Shannon Entropy) 
}
\usage{
clonality(x,entropy=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Numeric vector containing the counts of each species in a population}
  \item{entropy}{Logical indicating if the entropy should be reported, rather than the clonality (C=1-E)}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{cl}{Clonality}
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
