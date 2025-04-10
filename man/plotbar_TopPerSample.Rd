% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotbar_TopPerSample.R
\name{plotbar_TopPerSample}
\alias{plotbar_TopPerSample}
\title{plotbar_TopPerSample}
\usage{
plotbar_TopPerSample(
  ps,
  top = 10,
  relatief = TRUE,
  taxrank = "Species",
  taxfill = "Genus",
  output = "graph",
  x = "x_names",
  legend.position = "right",
  statistics = FALSE,
  GS = TRUE,
  angle = 0
)
}
\arguments{
\item{ps}{phyloseq object to be used}

\item{top}{the amount of most abundant ranks to be added in the plot}

\item{relatief}{Is the graph y-axis with realtive or absolute abundance (TRUE or FALSE)}

\item{taxrank}{the rank for the top. Mostly Species is used (is default)}

\item{taxfill}{the color in the graph}

\item{x}{the x-axis names. If unspecified the ps object names will be used}

\item{GS}{Make a GenusSpecies rank in the ps? (TRUE or FALSE). Default to TRUE, makes a GenusSpecies rank}

\item{angle}{angle for x-axis label (most likely 0 or 90°)}

\item{legend.position:}{Where is the legend? (left, right, bottom, top, or none). Legends can be very large}

\item{choosing}{"none" is equal to set legend=FALSE}

\item{statistics:}{Do yoiu want (relative) reads to be plotted (TRUE/FALSE)}

\item{output:}{undocumented feature what to export graph or table the graph is made of}
}
\description{
This Function plots a bar with a top(i) is per sample instead of an overall overall top(i).
Idea is to have higher explained reads compared to the ps_plotbar which uses overall top abundance. The functions returns a GGPLOT object
}
\examples{
plotbar_TopPerSample(ps, taxfill="Family", taxrank="Genus", top=10)
}
\keyword{plotbar}
