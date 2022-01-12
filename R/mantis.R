#' mantis.ScenarioOptions Returns a structure with all available Mantis Scenario options
#'
#'
#'
#' @return A data frame with the budget terms corresponding to the entire area
#' @export
#'
#' @examples
#' scenario<-mantis.ScenarioOptions()
mantis.ScenarioOptions <- function(){
  opt <-list(infile = 'incomingMSG.dat',
              outfile = 'testClientResults.dat',
              descr = 'This is a description of the simulation input',
              endSimYear = 2100,
              startRed = 2020,
              endRed = 2030,
              flowScen = 'C2VSIM_II_01',
              unsatScen = 'C2VSIM_SPRING_2000',
              wellType = 'VI',
              unsatWC = 0.01,
              rchScen = 'C2VSIM_II_01',
              bMap = 'Subregions',
              Regions = c('Subregion1','Subregion2'),
              RadSelect = c(), # X Y radius
              RectSelect = c(), # Xmin Ymin Xmax Ymax
              DepthRange = c(), #min max
              ScreenLenRange = c(), #min max
              loadScen = 'GNLM',
              loadSubScen = '',
              modifierName = '',
              modifierType = '',
              isLoadConc = 1,
              Crops = matrix(c(-9,1),nrow = 1, ncol = 2),
              SourceArea = c(), #Npixels minPixels maxPixels percPixels
              PixelRadius = c(),
              DebugID = '')
    return(opt)
}



mantis.Run <- function(scenario, client_prog){
  # Remove the output file
  if (file.exists(scenario$outfile)){
    file.remove(scenario$outfile)
  }
  mantis.WriteInput(scenario)
}



mantis.WriteInput <-function(scenario){
  con <- file(scenario$infile, open = "w")
  write(paste("endSimYear", scenario$endSimYear), file = con)
  write(paste("startRed", scenario$startRed), file = con, append = T)
  write(paste("endRed", scenario$endRed), file = con, append = T)
  write(paste("flowScen", scenario$flowScen), file = con, append = T)
  write(paste("unsatScen", scenario$unsatScen), file = con, append = T)
  write(paste("wellType", scenario$wellType), file = con, append = T)
  write(paste("unsatWC", scenario$unsatWC), file = con, append = T)
  write(paste("rchScen", scenario$rchScen), file = con, append = T)
  write(paste("bMap", scenario$bMap), file = con, append = T)

}
