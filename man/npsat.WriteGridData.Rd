% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npsat.R
\name{npsat.WriteGridData}
\alias{npsat.WriteGridData}
\title{npsat.WriteGridData writes data values for gridded interpolation}
\usage{
npsat.WriteGridData(filename, data, method, axisFiles)
}
\arguments{
\item{filename}{This is the name of the file}

\item{data}{The data to print. this is either an 1D vector, a 2D matrix or 3D array}

\item{method}{METHOD is LINEAR or NEAREST}

\item{axisFiles}{is the names of the files that describe the axis objects.}
}
\value{
just prints a file
}
\description{
npsat.WriteGridData writes data values for gridded interpolation
}
\examples{
2D Example:
N = 20
P <- peaks(v = N)
x <- 10*P$X[1, ]
y <- 10*P$Y[, 1]
orig <- -30
dx <- diff(x)[1]
writeAxis("peakAxis_cnst.tmp", c(orig, dx, N), "CONST")
writeData("peak_data_cnst_p.tmp", data = P$Z, mode = "POINT", "LINEAR",
           axisFiles = c("Rgridinterp/peakAxis_cnst.tmp", "Rgridinterp/peakAxis_cnst.tmp"))
Note that we pass as many axisFiles as needed. In this example the x and y axis are identical so we pass the same file
}
