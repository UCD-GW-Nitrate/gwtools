% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npsat.R
\name{npsat.WriteScattered}
\alias{npsat.WriteScattered}
\title{npsat.WriteScattered writes to a file the}
\usage{
npsat.WriteScattered(filename, PDIM, TYPE, MODE, DATA)
}
\arguments{
\item{filename}{The name of the file}

\item{PDIM}{is the number of columns that correspond to coordinates.  This is 1 for 1D points of 2 for 2D points.}

\item{TYPE}{Valid options for type are FULL, HOR or VERT}

\item{MODE}{Valid options for mode are SIMPLE or STRATIFIED}

\item{DATA}{the data to be printed. Data should have as many columns as needed.
For example it can be :
[x v]
[x v1 z1 v2 z2 ... vn-1 zn-1 vn]
[x y v]
[x y v1 z1 v2 z2 ... vn-1 zn-1 vn]
For 2D interpolation such as top, bottom elevation or recharge
npsat.WriteScattered(filename, 2, "HOR", "SIMPLE", data)}
}
\description{
npsat.WriteScattered writes to a file the
}
