\name{morisita}
\alias{morisita}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the Morisita's Overlap between two samples
}
\description{
The Morisita Overlap Index is a weighted overlap metric which takes into account the relative abundance of each species, as well as the size and Simpson's Diversity of each sample.  
}
\usage{
morisita(x,y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Numeric vector containing the counts of each species in sample 1}
  \item{y}{Numeric vector containing the counts of each species in sample 2}
}
\details{
x and y must have the same length.
}
\value{
\item{m}{Morista's Overlap Index}
}
\references{
https://en.wikipedia.org/wiki/Morisita's_overlap_index
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
