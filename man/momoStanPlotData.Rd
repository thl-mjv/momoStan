% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools.R
\name{momoStanPlotData}
\alias{momoStanPlotData}
\alias{momoStanFeature}
\alias{momoStanPlot}
\title{Plot methods for momoStanEffs}
\usage{
momoStanPlotData(obj, type = "baseline", colour = "red", group = 1,
  alpha = 0.2)

momoStanFeature(plt, obj, type = "baseline", colour = "red", group = 1,
  alpha = 0.2)

momoStanPlot(obj, type = "baseline", colour = "red", group = 1,
  alpha = 0.2)
}
\arguments{
\item{obj}{an object to plot}

\item{type}{type of the plot}

\item{colour}{colour of the ribbon}

\item{group}{which group to be plotted}

\item{alpha}{transparency of the ribbon}

\item{plt}{a plot to be added to}
}
\value{
a data or a ggplot
}
\description{
Plot methods for momoStanEffs
}
\section{Functions}{
\itemize{
\item \code{momoStanFeature}: additional features to a plot

\item \code{momoStanPlot}: plot data from momoStanEffs
}}

