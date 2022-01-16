
#' Convert subscripts to linear indices in column major order
#'
#' This function is somewhat the equivalent of the Matlab sub2ind.
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
#' This function is somewhat the equivalent of the Matlab ind2sub.
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


#' distPointLineSeg calculates the distance between a point and a line segment
#' defined by two points
#'
#' The function projects the point onto the line and returns the perpendicular distance
#' between the point and the line if the projection fall within the line segment.
#' If the projection falls outside the line segments it returns the closest distance between
#' the point and the closest end segment point
#'
#'
#' @param xp the x coordinate of the point
#' @param yp the y coordinate of the point
#' @param x1 the x coordinate of the end point of the line
#' @param y1 the y coordinate of the end point of the line
#' @param x2 the x coordinate of the other end point of the line
#' @param y2 the y coordinate of the other end point of the line
#'
#' @return the closest distance between the point and the line segment
#' @export
#'
#' @examples
distPointLineSeg <- function(xp, yp, x1, y1, x2, y2){
  L <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
  r <- ((y1 - yp) * (y1 - y2)-(x1 - xp)*(x2 - x1))/L^2
  if (r > 1 | r < 0){
    dst = min(sqrt( (xp - x1)^2 + (yp - y1)^2 ), sqrt( (xp - x2)^2 + (yp - y2)^2 ))
  }
  else{
    dst <- abs(((x2 - x1) * (y1 - yp) - (x1 - xp) * (y2 - y1))/ L)
  }
  return(dst)
}
