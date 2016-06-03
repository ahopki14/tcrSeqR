\name{iseqr_plot_factor}
\alias{iseqr_plot_factor}
\title{Convenience function for a plot_ds object}
\description{Simplifies plotting an immunoSeqR plot_ds object. The function returns a 
dotplot of the data in plot_ds$<metric>, splitting the data by whatever levels are in 
plot_ds$<x_val>, and optionally restricting the data to plot_ds$type==<type>.}

\usage{g <- iseqr_plot_factor(plot_ds,metric,x_val,type)}
\arguments{
  \item{plot_ds}{An immunoSeqR plot_ds object}
  \item{metric}{Which metric to be plotted}
  \item{x_val}{The metadata field on which to split the data for plotting. This must be a 
factor.}
  \item{type}{An optional argument to restrict the plot to a certain value in the type 
column}
}
\details{}
\value{
  \item{g}{A ggplot object containing the plot}
}
\references{}
\author{
Alexander Hopkins
}
\note{}

\seealso{iseqr_plot_metric}
\examples{}
\keyword{}
