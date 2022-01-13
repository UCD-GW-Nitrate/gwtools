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
              bMap = 'Subregions',
              Regions = c('Subregion1','Subregion2'),
              RadSelect = NA, # X Y radius
              RectSelect = NA, # Xmin Ymin Xmax Ymax
              DepthRange = NA, #min max
              ScreenLenRange = NA, #min max
              loadScen = 'GNLM',
              loadSubScen = NA,
              modifierName = NA,
              modifierType = NA,
              modifierUnit = NA,
              Crops = matrix(c(-9,1),nrow = 1, ncol = 2),
              SourceArea = NA, #Npixels minPixels maxPixels percPixels
              PixelRadius = NA,
              DebugID = NA)
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
  scenNames <- names(scenario)
  appnd <- FALSE
  if (!is.na(scenario$descr)){
    write(paste('#', scenario$descr) , file = con, append = appnd)
    appnd <- TRUE
  }

  for (i in 4:length(scenNames)) {
    if (!is.na(scenario[[i]][1])){
      if (scenNames[i] == 'Regions'){
        tmp <- paste('Nregions', length(scenario$Regions))
        write(paste(tmp, paste(scenario$Regions, collapse= " ") ) , file = con, append = appnd)
      }
      else if (scenNames[i] == 'Crops'){
        write(paste('Ncrops', dim(scenario$Crops)[1]) , file = con, append = appnd)
        for (j in 1:dim(scenario$Crops)[1]) {
          write(paste(scenario$Crops[j,1], scenario$Crops[j,2]) , file = con, append = appnd)
        }
      }
      else{
        tmp <- paste(scenario[[i]],collapse= " ")
         write(paste(scenNames[i], tmp) , file = con, append = appnd)
      }
      if (!appnd){
        appnd <- TRUE
      }
    }
  }
  close(con)
}

mantis.ReadOutput <- function(filename){
  con <- file(filename, open = "r")
  line <- readLines(con)
  close(con)
  tmp <- utils::read.table(file='filename',skip = 1)
}
