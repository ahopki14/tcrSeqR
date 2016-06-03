\name{find_clones}
\alias{find_clones}
\title{Finds clones of interest in an immunoSeqR dataset}
\description{Reads a list of clones from a text file, finds them in the dataset and then 
returns the restricted dataset.}
\usage{ds_coi <- find_clones(path,ds)}
\arguments{
  \item{path}{Path to a text file with one amino acid sequence per line}
  \item{ds}{An immunoSeqR dataset}
}
\details{}
\value{
  \item{ds_coi}{A subset of the input dataset, containing only the clones of interest}
}
\references{}
\author{
Alexander Hopkins
}
\note{}

\seealso{}
\examples{}
\keyword{}
