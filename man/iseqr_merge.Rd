\name{iseqr_merge}
\alias{iseqr_merge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Merge exported ImmunoSeq data files into a data frame
}
\description{

}
\usage{
iseqr_merge(all_files)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{all_files}{Character vector of paths to files to be merged. Should also provide unique information about each samples. This list will be sanitized and used as the sample names in the data.frame}
}
\details{

}
\value{
\item{ds}{data.frame containing counts. Column headers are replaced with names derived from all_files.}
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
\examples{}

