#' iwfm.readGWZBUD Reads the groundwater zone budget hdf5 file
#'
#' @param filename is the name of the file.
#' @param comp is a keyword of the groundwater component to extract.
#'   DP: Deep Percolation
#'   BS: Beginning Storage
#'   ES: Ending Storage
#'   NDP: Net Deep Percolation
#'   GFS: Gain From Stream
#'   R: Recharge
#'   GFL: Gain From Lake
#'   BI: Boundary Inflow
#'   S: Subsidence
#'   SI: Subsurface Irrigation
#'   TDO: Tile Drain Outflow
#'   P: Pumping
#'   NSI: Net Subsurface Inflow
#'   D: Discrepancy
#'   CS: Cumulative Subsidence
#'
#' @return A list with the groundwater budget components
#' @export
#'
#' @examples
#' GWZBUD <- iwfm.readGWZBUD("C2VSimFG_GW_ZBudget.hdf")
iwfm.readGWZBUD <- function(filename, comp = "ES"){

  return(NULL)
}
