#' nsgaii.readParetoHistory
#' Reads the pareto history file
#'
#' @param filename The file where the Pareto history is written
#' @param nGen Then total number of generations
#'
#' @return a list of lenght nGew with a matrix in each entry representing the
#'   current pareto front
#' @export
nsgaii.readParetoHistory <- function(filename, nGen){
  l <- vector(mode = "list", length = nGen)
  nlinesskip <- 0
  for (i in 1:nGen) {
    info <- scan(file = filename, skip = nlinesskip, n = 3)
    nlinesskip = nlinesskip + 1
    p <- scan(file = filename, skip = nlinesskip, n = info[2]*info[3])
    nlinesskip = nlinesskip + info[3]
    m <- matrix(data = p, nrow = info[3], ncol = info[2], byrow = TRUE)
    l[[i]] <- m
  }
  return(l)
}


#' nsgaii.readParetoSolution reads the pareto solution.
#'
#'
#' @param filename The file name is the OutputFile parameter of the configuration file of nsgaii
#'
#' @return a list with 2 elements. dv contains the decision variables and of the
#' objective function values
#'
#' @export
nsgaii.readParetoSolution <- function(filename){
  l <- vector(mode = "list", length = 2)
  info <- scan(file =  filename, n = 3)
  d <- scan(file = filename, n = info[1]*info[2], skip = 1)
  of <- scan(file = filename, n = info[1]*info[3],skip = info[1]+1)
  dvmat <- matrix(data = d, nrow = info[1], ncol = info[2], byrow = TRUE)
  ofmat <- matrix(data = of, nrow = info[1], ncol = info[3], byrow = TRUE)
  return(list(dv = dvmat, of = ofmat))
}

#' nsgaii.calculateHyperVolume Calculates the hypervolume of a given pareto front
#' with respect to a reference point
#'
#' @param pareto This is the pareto front
#' @param RefPnt This is the reference point. That is a point that is dominated by every solution
#' in the pareto front
#' @export
nsgaii.calculateHyperVolume <- function(pareto, RefPnt){
  psort <- sort(pareto[,1], index.return = TRUE)

  ob1 <- psort$x
  ob2 <- pareto[psort$ix,2]
  vol <- abs(ob1[1] - RefPnt[1]) * abs(ob2[1] - RefPnt[2])
  for (i in 2:length(ob1)) {
    vol <- vol + abs(ob1[i] - RefPnt[1])*abs(ob2[i] - ob2[i-1])
  }
  return(vol)
}
