# How to build and install the package
# devtools::document()
# devtools::load_all()
# devtools::install()

#' Writes a mesh into a file for reading from the npsat_engine
#'
#'
#' @param filename is the name of the file to write the mesh.
#' @param nd A matrix n x 2 with the mesh node coordinates
#' @param msh A matrix with the indices of the mesh elements. The matrix must be based on 1.
#' However the file will be printed as zero based
#' @return nothing
#' @export
npsat.writeMesh <- function(filename, nd, msh){
  write(c(dim(nd)[1], dim(msh)[1]), file = filename, append = FALSE)
  utils::write.table(nd, file = filename, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
  utils::write.table(msh, file = filename, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
}


#' npsat.WriteScattered writes to a file the
#'
#' @param filename The name of the file
#' @param PDIM is the number of columns that correspond to coordinates.  This is 1 for 1D points of 2 for 2D points.
#' @param TYPE Valid options for type are FULL, HOR or VERT
#' @param MODE Valid options for mode are SIMPLE or STRATIFIED
#' @param DATA the data to be printed. Data should have as many columns as needed.
#' For example it can be :
#' [x v]
#' [x v1 z1 v2 z2 ... vn-1 zn-1 vn]
#' [x y v]
#' [x y v1 z1 v2 z2 ... vn-1 zn-1 vn]
#'
#' @return
#' @examples
#' For 2D interpolation such as top, bottom elevation or recharge
#' npsat.WriteScattered(filename, 2, "HOR", "SIMPLE", data)
#'
#' @export
npsat.WriteScattered <- function(filename, PDIM, TYPE, MODE, DATA){
  write("SCATTERED", file = filename, append = FALSE)
  write(TYPE, file = filename, append = TRUE)
  write(MODE, file = filename, append = TRUE)
  Ndata <- dim(DATA)[2] - PDIM
  write(paste(dim(DATA)[1], Ndata), file = filename, append = TRUE)
  utils::write.table(DATA, file = filename, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
}

#' npsat.ReadScattered reads the scattered format files
#'
#' @param filename is the name of the file
#' @param PDIM is the dimention of the point. Valid options are 2 for 2D points or 1 for 1D points.
#' At the moment the points can be only 1D or 2D.
#'
#' @return a data frame with the X Y Values
#' @export
npsat.ReadScattered <- function(filename, PDIM = 2){
  # Read header
  tmp <- readLines(filename, n = 4)
  N <- as.numeric(strsplit(tmp[4], " ")[[1]])
  if (PDIM == 1){
    cnames <- c("X")
  }
  else if (PDIM == 2){
    cnames <- c("X", "Y")
  }
  if(N[2] == 1){
    cnames <- c(cnames, "V")
  }
  else if (N[2] > 1){
    cnames <- c(cnames, "V1")
    for (i in seq(1, (N[2]-1)/2, 1)) {
      cnames <- c(cnames, paste0("Z", i), paste0("V", i+1))
    }

  }
  DATA <- utils::read.table(file = filename,
                     header = FALSE, sep = "", skip = 4, nrows = N[1],
                     quote = "",fill = TRUE,
                     col.names = cnames)

  return(DATA)
}

#' npsat.WriteWells writes the wells to a file
#'
#' @param filename is the name of the file
#' @param wells is a data frame with the following files in the following order
#' X coordinate
#' Y coordinate
#' Top elevation of the top of the screen
#' Bottom elevation of the bottom of the screen
#' Q pumping rate. Extraction is negative
#'
#' @return nothing
#' @export
npsat.WriteWells <- function(filename, wells){
  write(dim(wells)[1], file = filename, append = FALSE)
  utils::write.table(wells, file = filename, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
}


#' npsat.ReadStreams reads the Stream file
#'
#' @param filename is the name of the file
#' @param maxChar is the maximum number of characters in a single line
#'
#' @return A list of Stream polygons with the following fields:
#' Coords: the coordinates of the stream polygon
#' Rate: The rate in L/T units
#' Area the area of the polygon
#' Q: The total volume of water in L^3/T
#' Info: This can be any info associated with the river segments
#'
#' @export
npsat.ReadStreams <- function(filename, maxChar = 1000){
  allLines <- readLines(filename)
  # read the number of polygons
  ind <- 1
  tmp <- strsplit(substr(allLines[ind], 1, maxChar)[[1]], split = " ")[[1]]
  tmp <- tmp[which(tmp != "")]
  Npoly <- as.numeric(tmp)
  streamCoords <- vector(mode = "list", length = Npoly)
  Area <- vector(mode = "numeric", length = Npoly)
  Rate <- vector(mode = "numeric", length = Npoly)
  Q <- vector(mode = "numeric", length = Npoly)
  Info <- vector(mode = "list", length = Npoly)
  for (i in 1:Npoly) {
    ind <- ind + 1
    tmp <- strsplit(substr(allLines[ind], 1, maxChar)[[1]], split = " ")[[1]]
    tmp <- tmp[which(tmp != "")]
    tmp <- as.numeric(tmp)
    Npnts <- tmp[1]
    Rate[i] <- tmp[2]
    if (length(tmp) > 2){
      Info[[i]] <- tmp[-c(1,2)]
    }
    coords <- matrix(data = NA, nrow = Npnts, ncol = 2)
    for (j in 1:Npnts) {
      ind <- ind + 1
      tmp <- strsplit(substr(allLines[ind], 1, maxChar)[[1]], split = " ")[[1]]
      tmp <- tmp[which(tmp != "")]
      tmp <- as.numeric(tmp)
      coords[j,] <- tmp
    }
    Area[i] <- pracma::polyarea(coords[,1], coords[,2])
    Q[i] <- Area[i]*Rate[i]
    streamCoords[[i]] <- coords
  }
  return(list(Coords = streamCoords, Rate = Rate, Area = Area, Q = Q, Info = Info))
}

#' npsat.ReadWells reads the well file
#'
#' @param filename is the name of the file
#' @param maxChar is the maximum number of characters in a single line
#'
#' @return A dataframe with the wells with the following fields:
#' \itemize{
#'     \item X, Y: coordinates of the well
#'     \item T: top elevation of the well screen
#'     \item B: bottom elevation of the well screen
#'     \item Q: The pumping rate in L^3/T
#' }
#'
#' @export
npsat.ReadWells <- function(filename){
  M <- utils::read.table(file = filename,
             header = FALSE, sep = "", skip = 1,
             quote = "",fill = TRUE,
             col.names = c("X", "Y", "T", "B", "Q"))
  return(M)
}


#' npsat.ReadWaterTableCloudPoints reads the water table cloud point file
#'
#' @param prefix is the prefix used in the simulation input file
#' @param suffix This is 'top' by default
#' @param nproc The number of processors used in the simulation
#' @param iter the iteration number
#' @print_file This is currently not used
#'
#' @return A dataframe with the following fields:
#' \itemize{
#'     \item X, Y: coordinates of the node
#'     \item Zold: the initial elevation
#'     \item Znew: the new elevation
#' }
#'
#' @export
npsat.ReadWaterTableCloudPoints <- function(prefix, suffix = 'top', nproc = 1, iter = 1, print_file = F){
  iter <- iter - 1
  df <- data.frame(matrix(data = NA, nrow = 0, ncol = 4))
  for (i in 0:(nproc-1)) {
    fname <- paste0(prefix,suffix,'_',sprintf("%03d",iter),'_', sprintf("%04d",i), '.xyz')
    npoints <- read.table(file = fname,
                          header = FALSE, sep = "", skip = 0, nrows = 1,
                          quote = "",fill = TRUE,
                          col.names = c("N"))
    if (npoints == 0){
      next
    }
    datapnts <- read.table(file = fname,
                           header = FALSE, sep = "", skip = 1, nrows = npoints$N,
                           quote = "",fill = TRUE)
    df <- rbind(df, datapnts)
  }
  x <- c("X", "Y", "Zold", "Znew")
  colnames(df) <- x
  return(df)
}


#' npsat.WriteGridAxis writes a xis object to file
#'
#' @param filename is the file name
#' @param x If the spacing option is CONST then x is a length 3 vector
#' where the x[1] is the origin x[2] is the distance between ticks and x[3] is the number of ticks
#' in that axis. If spacing is VAR then x must be a monotonous increasing vector with the ticks.
#' This is a good option if the ticks are not evenly spaced
#' @param spacing spacing is either CONST or VAR
#'
#' @return just prints a file
#' @export
#'
#' @examples
#' npsat.WriteGridAxis("axis1.dat, c(0,100,23), "CONST")
#' will generate an axis with ticks starting from zero to 2200
#' npsat.WriteGridAxis("axis1.dat, x, "VAR")
#' will print the vector x as the axis. Note that it is not checked for doubles or
#' non monotonous increase in this function neither by the program itself
npsat.WriteGridAxis <- function(filename, x, spacing){
  if (spacing == "CONST"){
    N <- x[3]
  }
  else{
    N <- length(x)
  }
  write(paste(spacing,N), file = filename, append = F)
  if (spacing == "CONST"){
    write(x[1:2], file = filename, append = T)
  }
  else{
    write.table(x, file = filename, append = T, row.names = F, col.names = F)
  }
}

#' npsat.WriteGridData writes data values for gridded interpolation
#'
#' @param filename This is the name of the file
#' @param data The data to print. this is either an 1D vector, a 2D matrix or 3D array
#' @param method METHOD is LINEAR or NEAREST
#' @param axisFiles is the names of the files that describe the axis objects.
#'
#' @return just prints a file
#' @export
#'
#' @examples
#' 2D Example:
#' N = 20
#' P <- peaks(v = N)
#' x <- 10*P$X[1, ]
#' y <- 10*P$Y[, 1]
#' orig <- -30
#' dx <- diff(x)[1]
#' writeAxis("peakAxis_cnst.tmp", c(orig, dx, N), "CONST")
#' writeData("peak_data_cnst_p.tmp", data = P$Z, mode = "POINT", "LINEAR",
#'            axisFiles = c("Rgridinterp/peakAxis_cnst.tmp", "Rgridinterp/peakAxis_cnst.tmp"))
#' Note that we pass as many axisFiles as needed. In this example the x and y axis are identical so we pass the same file
npsat.WriteGridData <- function(filename, data, method, axisFiles){
  if (is.null(dim(data))){
    write(paste(method, length(data)), file = filename, append = F)
    write(axisFiles, file = filename, append = T, ncolumns = length(axisFiles))
    write.table(data, file = filename, append = T, row.names = F, col.names = F)
  }
  else{
    if (length(dim(data)) == 2){
      write(paste(method, dim(data)[2], dim(data)[1]), file = filename, append = F)
      write(axisFiles, file = filename, append = T, ncolumns = length(axisFiles))
      write.table(data, file = filename, append = T, row.names = F, col.names = F)
    }
    else if (length(dim(data)) == 3){
      write(paste(method, dim(data)[2], dim(data)[1], dim(data)[3]), file = filename, append = F)
      write(axisFiles, file = filename, append = T, ncolumns = length(axisFiles))
      for (i in 1:dim(data)[3]) {
        write.table(data[,,i], file = filename, append = T, row.names = F, col.names = F)
      }
    }
  }
}

npsat.WriteWells4ichnos <- function(filename, wells){

}
