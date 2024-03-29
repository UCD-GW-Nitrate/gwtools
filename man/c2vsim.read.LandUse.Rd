% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/c2vsim.R
\name{c2vsim.read.LandUse}
\alias{c2vsim.read.LandUse}
\title{c2vsim.read.LandUse
Reads the land use time series. For the IWFM version 3 it reads the CVland use.
For the IWFM 15 it reads the files that have the _area suffix.
The default values correspond to the coarse grid case version 3}
\usage{
c2vsim.read.LandUse(
  filename,
  NtimeSteps = 88,
  Nelem = 1392,
  Ninfo = 5,
  colNames = NA,
  Nskip = 95,
  maxChar = 1000
)
}
\arguments{
\item{filename}{The name of the file to read}

\item{NtimeSteps}{The number of time steps in the file}

\item{Nelem}{The number of elements}

\item{Ninfo}{The number of columns to read}

\item{Nskip}{The number of lines to skip before the timeseries}

\item{maxChar}{The number of characters to read for each line}
}
\value{
The Land use
}
\description{
c2vsim.read.LandUse
Reads the land use time series. For the IWFM version 3 it reads the CVland use.
For the IWFM 15 it reads the files that have the _area suffix.
The default values correspond to the coarse grid case version 3
}
