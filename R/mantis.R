# How to build and install the package
# devtools::document()
# devtools::load_all()
# devtools::install()

#' mantis.ScenarioOptions Returns a structure with all available Mantis Scenario options
#'
#' @return A data frame with the budget terms corresponding to the entire area
#' @export
#'
#' @examples
#' scenario<-mantis.ScenarioOptions()
mantis.ScenarioOptions <- function(){
  opt <-list( client = NA,
              infile = 'incomingMSG.dat',
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
              minRch = 10,
              maxConc = -1,
              loadTrans = NA,
              loadTransYearStart = NA,
              loadTransYearEnd = NA,
              Crops = matrix(c(-9,1),nrow = 1, ncol = 2),
              SourceArea = NA, #Npixels minPixels maxPixels percPixels
              PixelRadius = NA,
              DebugID = NA)
    return(opt)
}


#' mantis.Run Runs the scenario and returns the results
#'
#' @param scenario is a list with the scenario options
#'
#' @return a list with the well breakthrough curves
#' @export
#'
#' @examples
#' scenario <- mantis.ScenarioOptions()
#' scenario$client <- 'MantisClient.exe'
#' result <- mantis.Run(scenario)
mantis.Run <- function(scenario){
  if (is.na(scenario$client)){
    outmsg <- 'Client Path not specified'
    return(outmsg)
  }

  # Remove the output file
  if (file.exists(scenario$outfile)){
    file.remove(scenario$outfile)
  }
  mantis.WriteInput(scenario)

  system2(command = scenario$client, args = c(scenario$infile, scenario$outfile))

  res <- mantis.ReadOutput(scenario$outfile)
  return(res)
}


#' mantis.Quit Stops the server
#'
#' @param scenario the only useful option for this function
#' is the client which is is the path to the client executable
#'
#' @return Returns nothing. It just quits the server
#' @export
#'
#' @examples
#' scenario <- mantis.ScenarioOptions()
#' scenario$client <- 'MantisClient.exe'
#' mantis.Quit <- function( scenario)
mantis.Quit <- function(scenario){
  if (!is.na(scenario$client)){
    system2(command = scenario$client, args = c('quit'))
    return(1)
  }
  else{
    outmsg <- 'Client Path not specified'
    return(outmsg)
  }
}

#' mantis.Ping Test the communication with server.
#' It sends a ping and expects a pong
#'
#' @param scenario the only useful option for this function
#' is the client which is is the path to the client executable
#'
#' @return Returns a pong if the communication is succesful
#' @export
#'
#' @examples
#' scenario <- mantis.ScenarioOptions()
#' scenario$client <- 'MantisClient.exe'
#' mantis.Ping <- function( scenario)
mantis.Ping <- function(scenario){
  if (!is.na(scenario$client)){
    system2(command = scenario$client, args = c('ping', scenario$outfile))
    outmsg <- mantis.ReadOutput(scenario$outfile)
  }
  else{
    outmsg <- 'Client Path not specified'
  }
  return(outmsg)
}


#' mantis.WriteInput writes the input file for a given scenario.
#' This is called by the mantis.Run and there is no reason why
#' someone would need to run this
#'
#' @param scenario is a list with the scenario options
#' @param client_prog is the path to the Mantis client executable
#'
#' @return It returns nothing and prints the input file.
#' The name of the input file is set by the option scenario$infile
#' @export
#'
#' @examples
#' mantis.WriteInput <-function(scenario)
mantis.WriteInput <-function(scenario){
  con <- file(scenario$infile, open = "w")
  scenNames <- names(scenario)
  appnd <- FALSE
  if (!is.na(scenario$descr)){
    write(paste('#', scenario$descr) , file = con, append = appnd)
    appnd <- TRUE
  }

  for (i in 5:length(scenNames)) {
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


#' mantis.ReadOutput reads the output file of the mantis client.
#' This is called by the mantis.Run and there is no reason why
#' someone would need to run this
#'
#' @param filename This is the output file. The file name is specified
#' with the option scenario$outfile
#'
#' @return returns a table with the wells break though curves.
#' The number of rows is equal to the number of wells in the selected area
#' and the number of columns is the number of simulation years
#' @export
#'
#' @examples
#' mantis.ReadOutput <- function(scenario$outfile)
mantis.ReadOutput <- function(filename){
  con <- file(filename, open = "r")
  line <- readLines(con,1)
  close(con)
  t <- strsplit(substr(line,1,100), split = " ")
  NNs <- as.numeric(t[[1]])
  if (NNs[1] == 0){
    alllines <- readLines(filename)
    return(alllines)
  }
  else{
    res <- utils::read.table(file = filename, skip = 1)
    return(res)
  }
}
