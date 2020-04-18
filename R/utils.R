
#' Convert subscripts to linear indices in column major order
#'
#' This function is somewhat the equivelant of the matlab sub2ind.
#'
#' @param i row index
#' @param j column index
#' @param n number of rows
#' @return the linear index.
#' @usage ind <- sub2ind(i,j,n), then A[i,j] == A[ind] if dim(A)[1] == n
#' @export
sub2ind <- function(i,j,n){
  n*(j-1)+i
}

#' Convert linear indices to subscripts in column major order
#'
#' This function is somewhat the equivelant of the matlab ind2sub.
#'
#' @param l linear index
#' @param n number of rows
#' @return the linear index.
#' @usage t <- ind2sub(j,n), then A[t$r, t$c] == A[l] if dim(A)[1] == n
#' @export
#' @export
ind2sub <- function(l,n){
  c <- ceiling (l / n)
  r <- l - (c-1)*n
  return(cbind(r,c))
}
