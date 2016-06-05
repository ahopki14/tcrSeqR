\name{iseqr_lookup}
\alias{iseqr_lookup}
\title{A function to pull metadata from an immunoSeqR dictionary}
\description{This is a function to easily query a dictionary for a single piece of metadata.}
\usage{out <- iseqr_lookup(i='2',dict,i_col='patient',o_col='response'}
\arguments{
  \item{i}{The input of the query}
  \item{dict}{An immunoSeqR dictionary}
  \item{i_col}{The column in which to look for <i>}
  \item{o_col}{The output column}
}
\details{}
\value{
  \item{out}{The output of the query}
}
\references{}
\author{
Alexander Hopkins
}
\note{}

\seealso{}
\examples{
# retrieve the response status of patient 2
out <- iseqr_lookup(i='2',dict,i_col='patient',o_col='response')

# retrieve the treatment arm of patient 2
out <- iseqr_lookup(i='2',dict,i_col='patient',o_col='arm')

# make a vector of all the responses
patients <- levels(dict$patient)
out <- vector()
for(i in patients){
  out[i] <- iseqr_lookup(i,dict)
}

}
\keyword{}
