\name{iseqr_clean_ds}
\alias{iseqr_clean_ds}
\title{Clean Dataset After Subsetting}
\description{Removes any rows of a dataset with no counts in any sample. This is used after 
restricting an immunoSeqR dataset to a smaller subset of samples.}
\usage{ds <- iseqr_clean_ds(ds)}
\arguments{
  \item{ds}{An immunoSeqR dataset}
}
\details{}
\value{
  \item{ds}{A cleaned version (with potentially fewer rows) of the input dataset}
}
\references{}
\author{
Alexander Hopkins
}
\note{It does not refactor the aa or nt columns. }

\seealso{}
\examples{}
\keyword{}
