\name{iseqr_plot_metric}
\alias{iseqr_plot_metric}
\title{Convenience function for a plot_ds object}
\description{Simplifies plotting an immunoSeqR plot_ds object. The function returns a 
scatter of the data in plot_ds$<metric> vs. plot_ds$<x_val>, and optionally restricting the 
data to plot_ds$type==<type>.}

\usage{g <- iseqr_plot_metric(plot_ds,metric,x_val,type,sm)}
\arguments{
  \item{plot_ds}{An immunoSeqR plot_ds object}
  \item{metric}{Which metric to be plotted}
  \item{x_val}{The metadata field to plot the data against. This must be a numeric vector.}
  \item{type}{An optional argument to restrict the plot to a certain value in the type 
column}
  \item{sm}{Logical indicating if the geom_smooth linear model should also be plotted.}
  \item{hjust}{Distance to place text labels from point. Default is -0.5}
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

\seealso{iseqr_plot_factor}
\examples{}
\keyword{}
