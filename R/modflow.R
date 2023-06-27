#' modflow.readCBC reads the cell by cell modflow binary file
#' See more on https://water.usgs.gov/nrp/gwsoftware/modflow2000/MFDOC/index.html?frequently_asked_questions.htm
#' expand L and then expand Flow Data
#'
#' @param filename is the file name
#' @param endian should be either "little" or "big"
#' @param intsize number of bytes for the integer numbers
#'
#' @return A complicated List
#' @export
#'
#' @examples cbc <- modflow.readCB(cbc.out)
modflow.readCBC <- function(filename, endian = "little", intsize = 4){
  flbn <- file(filename,"rb")
  Out <- vector(mode = "list", length = 0)
  dataList <- vector(mode = "list", length = 0)
  TimeStepList <- vector(mode = "list", length = 4)
  KPERprevious = -9;

  while(T){
    tmp_data <- NA
    KSTP <- readBin(flbn, what = "int", endian = endian,  size = intsize)
    if (length(KSTP) == 0){break}

    KPER <- readBin(flbn, what = "int", endian = endian,  size = intsize)
    if (length(KPER) == 0){break}

    TMP <- readBin(flbn, what = "int", endian = endian, n=16, size = 1)
    if (length(TMP) == 0){break}
    DESC <- intToUtf8(TMP)
    DESC <- trimws(DESC)

    NCOL <- readBin(flbn, what = "int", endian = endian,  size = intsize)
    if (length(NCOL) == 0){break}
    NROW <- readBin(flbn, what = "int", endian = endian,  size = intsize)
    if (length(NROW) == 0){break}
    NLAY <- readBin(flbn, what = "int", endian = endian, size = intsize)
    if (length(NLAY) == 0){break}

    if (NLAY < 0){
      NLAY <- abs(NLAY)
      ITYPE <- readBin(flbn, what = "int", endian = endian, size = intsize)
      DELT <- readBin(flbn, what = "numeric", endian = endian, size = 4)
      PERTIM <- readBin(flbn, what = "numeric", endian = endian, size = 4)
      TOTIM <- readBin(flbn, what = "numeric", endian = endian, size = 4)
      if (ITYPE == 0 || ITYPE == 1){
        tmp_data <- readBin(flbn, what = "numeric", endian = endian, n = NLAY*NCOL*NROW, size = 4)
        tmp_data <- array(data = tmp_data, dim = c(NCOL, NROW, NLAY))
        tmp_data <- aperm(tmp_data, c(2,1,3))
      }
      else if(ITYPE == 2 || ITYPE == 5){
        if (ITYPE == 5){
          NVAL <- readBin(flbn, what = "int", endian = endian, size = intsize)
          if (NVAL > 1){
            stop("ITYPE = 5 and NVAL > 1 is not implemented yet. (This is a good time to contact us at gkourakos@ucdavis.edu)")
          }
        }

        NLIST <- readBin(flbn, what = "int", endian = endian, size = intsize)
        if (NLIST > 0){
          NRC <- NROW*NCOL
          tmp_data <- matrix(data = NA, nrow = NLIST, ncol = 4)
          for (i in 1:NLIST) {
            ICELL <- readBin(flbn, what = "int", endian = endian, size = intsize)
            val <- readBin(flbn, what = "numeric", endian = endian, size = 4)
            lay <- floor((ICELL-1)/NRC+1);
            row <- floor(((ICELL-(lay-1)*NRC)-1)/NCOL+1);
            col <- ICELL-(lay-1)*NRC-(row-1)*NCOL;
            tmp_data[i,] <- c(lay, row, col, val)
          }
        }
      }
      else if (ITYPE == 3){
        tmp_data <- matrix(data = NA, nrow = NROW*NCOL, ncol = 2)
        tmp_data[,1] <- readBin(flbn, what = "int", endian = endian, n = NROW*NCOL, size = intsize)
        tmp_data[,2] <- readBin(flbn, what = "numeric", endian = endian, n = NCOL*NROW, size = 4)
      }
      else if (ITYPE == 4){
        stop("ITYPE: 4 is not implemented yet. (This is a good time to contact us at gkourakos@ucdavis.edu)")
      }
    }
    else if (NLAY > 0){
      tmp_data <- readBin(flbn, what = "numeric", endian = endian, n = NLAY*NCOL*NROW, size = 4)
      tmp_data <- array(data = tmp_data, dim = c(NCOL, NROW, NLAY))
      tmp_data <- aperm(tmp_data, c(2,1,3))
    }

    TimeStepList <- list(KSTP = KSTP,
                         KPER = KPER,
                         DESC = DESC,
                         DATA = tmp_data)

    if (KPER != KPERprevious){
      if (KPERprevious == -9){
        KPERprevious <- KPER
        print(KPER)
        dataList <- c(dataList, list(TimeStepList))
      }
      else{
        Out <- c(Out, list(dataList))
        dataList <- vector(mode = "list", length = 0)
        dataList <- c(dataList, list(TimeStepList))
        print(KPER)
        flush.console()
        KPERprevious <- KPER
      }
    }
    else{
      dataList <- c(dataList, list(TimeStepList))
    }

  }# while loop
  Out <- c(Out, list(dataList))
  close(flbn)
  return(Out)
}


#' modflow.readArrayASCII reads a 3D array that is printed in ASCII.
#' This can be used to read the modflow heads when they are printed as ASCII
#'
#' @param filename is the file name
#' @param nlay The number of layers
#' @param maxchar A guess about the maximum width of the lines
#'
#' @return A list with one 3D array per time step
#' @export
#'
#' @examples gwHead <- modflow.readArrayASCII(headsout.txt, nlay = 13)
modflow.readArrayASCII <- function(filename, nlay, maxchar = 1000){
  out <- vector(mode = "list", length = 0)
  info <- vector(mode = "list", length = 0)
  alllines <- readLines(filename)
  iline <- 1
  while(T){
    for (k in 1:nlay) {
      # Read the first line
      t <- substr(alllines[iline], 1, maxchar)
      iline <- iline + 1

      t <- strsplit(t, split = " ")
      t <- t[[1]][which(t[[1]] != "")]
      if (k == 1){
        print(paste("reading timestep: ",t[2]))
      }

      nrow <- as.numeric(t[7])
      ncol <- as.numeric(t[6])
      ilay <- as.numeric(t[8])
      datalay <- array(data = NA, dim = c(nrow,ncol))

      if (iline == 2){
        dataArray <- array(data = NA, dim = c(nrow, ncol, nlay))
      }
      datalay[] <- NA
      for (i in 1:nrow) {
        col_ind_s <- 1
        while(T){
          t <- substr(alllines[iline], 1, maxchar)
          iline <-  iline + 1
          t <- strsplit(t, split = " ")
          t <- as.numeric(t[[1]][which(t[[1]] != "")])
          col_ind_e <- col_ind_s + length(t)-1
          datalay[i,col_ind_s:col_ind_e] <- t
          col_ind_s <- col_ind_e + 1
          if (col_ind_s > ncol){
            break
          }
        }
      }
      dataArray[,,ilay] <- datalay
    }

    out <- c(out, list(dataArray))
    dataArray[] <- NA
    if (iline > length(alllines)){
      break
    }
  }
  return(out)
}
