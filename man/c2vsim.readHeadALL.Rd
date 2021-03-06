% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/c2vsim.R
\name{c2vsim.readHeadALL}
\alias{c2vsim.readHeadALL}
\title{c2vsim.readHeadALL Reads thegroundwater head for all layers and all time steps}
\usage{
c2vsim.readHeadALL(
  filename,
  nNode = 30179,
  nLay = 4,
  NtimeSteps = 505,
  nSkip = 6,
  quiet = FALSE
)
}
\arguments{
\item{filename}{is the name of the file with the groundwater head values.}

\item{nNode}{is the total number of nodes in the simulation.}

\item{nLay}{is the Number of layers.}

\item{NtimeSteps}{is the number of monthly time steps.}

\item{nSkip}{is the number of lines to skip before the first time step values are printed}
}
\value{
a list of lists. The first item in the list is a vector with the node ids.
The second lement in the list is a data frame with the time stamps
The third length NtimeSteps where each element in the list is a matrix
nNode x nLay.
}
\description{
The default values correspond to the fine grid version.
}
