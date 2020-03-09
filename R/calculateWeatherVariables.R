#' @title calculateDewPointTemp.
#'
#' @description
#' \code{calculateDewPointTemp} function for calculating Air dew point temperature
#'
#' @details
#' function for calculating Air dew point temperature
#' Equation 3-11 in FAO 56
#'
#' @references see aWhere API Documentation (insert HTTP address)
#'
#' @param - actualVapourPressure: double
#' 
#' @return double
#'
#' @examples
#' NEED TO INSERT
#' 
#' @export

calculateDewPointTemp <- function(actualVapourPressure){
      
      tDew <- (116.91 + 237.3 * log(actualVapourPressure))/
            (16.78-log(actualVapourPressure))
      
      return(tDew)
}

#' @title calculateSlopeOfVapourPressureCurve.
#'
#' @description
#' \code{calculateSlopeOfVapourPressureCurve} function for calcuating the slope of the vapour pressure curve 
#' 
#' @details
#' function for calcuating the slope of the vapour pressure curve based on temp
#' temp in C
#' eqn 13 from FAO 56
#'
#' @references see aWhere API Documentation (insert HTTP address)
#'
#' @param - temp: In C - double
#' 
#' @return double
#'
#' @examples
#' NEED TO INSERT
#' 
#' @export

calculateSlopeOfVapourPressureCurve <- function(temp) {
      
      slope <- (4098*(.6108*exp((17.27*temp)/(temp + 237.3))))/
                        (temp+237.3)^2
      
      return(slope)
      
}

#' @title calculateLatentHeatVap.
#'
#' @description
#' \code{calculateLatentHeatVap} function for calculating the latent heat of vaporization
#' 
#' @details
#' function for calculating the latent heat of vaporization
#' FAO 56 simplifies by using a single value over all temps as it changes little
#' 2.45 MJ kg-1.
#'
#' @references see aWhere API Documentation (insert HTTP address)
#' 
#' @return double
#'
#' @examples
#' NEED TO INSERT
#' 
#' @export

calculateLatentHeatVap <- function() {
      
      return(2.45)
}

#' @title calculateAtmosphericPressure.
#'
#' @description
#' \code{calculateAtmosphericPressure} #function for calculating Atmospheric Pressure
#' 
#' @details
#' function for calculating Atmospheric Pressure: equantion (3-2) in FAO 56
#' 
#' P atmospheric pressure at elevation z [kPa]
#' Po atmospheric pressure at sea level = 101.3 [kPa]
#' z elevation [m]
#' zo elevation at reference level [m] (normally 2)
#' g gravitational acceleration = 9.807 [m s-2]
#' R specific gas constant = 287 [J kg-1 K-1]
#' al constant lapse rate moist air = 0.0065 [K m-1]
#' TKo reference temperature [K] at elevation zo given by
#' TKo = 273.16 + T (3-3)
#' where: T mean air temperature for the time period of calculation [oC]
#'
#' @references see aWhere API Documentation (insert HTTP address)
#'
#' @param - elevation: In m - double
#' @param - elevAtRefLevel: In m - double
#' @param - temp: In C - double
#' 
#' @return double
#'
#' @examples
#' NEED TO INSERT
#' 
#' @export


calculateAtmosphericPressure <- function(elevation,elevAtRefLevel,temp){
      
      refTemp                       <- 273.16 + temp
      atmosphericPressureAtSeaLevel <- 101.3
      gravAcc                       <- 9.807
      specGasConstant               <- 287
      constLapseRateMoistAir        <- .0065
      
      atmosphericPressure <- 
            atmosphericPressureAtSeaLevel * ((refTemp - constLapseRateMoistAir * 
                      (elevation - elevAtRefLevel))/refTemp) ^ (gravAcc/(constLapseRateMoistAir * specGasConstant))
      
      return(atmosphericPressure)
      
}

#' @title calculateAtmosphericDensity.
#'
#' @description
#' \code{calculateAtmosphericDensity} function for calcuating the atmospheric density 
#' @details
#' function to calculate the atmospheric density from temp and atmospheric pressure
#' FAO 56 3-5 and 3-6
#'
#' @references see aWhere API Documentation (insert HTTP address)
#'
#' @param - temp: In C - double
#' @param - RH: In C - double
#' @param - elev: In C - double
#' 
#' @return double
#'
#' @examples
#' NEED TO INSERT
#' 
#' @export

calculateAtmosphericDensity <- function(temp,RH,elev) {
      
      specGasConstant     <- 287
      atmosphericPressure <- calculateAtmosphericPressure(elev,2,temp)
      actualVapPres       <- calculateActualVapourPressure(RH,RH,temp,temp)
      
      tempK            <- 273.16 + temp
      tempKV           <- tempK * (1-.378 * (actualVapPres/atmosphericPressure))^(-1)
      
      atmosphericDen <- (1000 * atmosphericPressure) / (tempKV * specGasConstant)
      
      return(atmosphericDen)
}


#' @title calculateSaturationVapourPressure.
#'
#' @description
#' \code{calculateSaturationVapourPressure} function for calcuating the saturation vapour pressure
#' 
#' @details
#' function to calculate the saturation vapour pressure
#' FAO 3-8, temp in C
#'
#' @references see aWhere API Documentation (insert HTTP address)
#'
#' @param - temp: In C - double
#' 
#' @return double
#'
#' @examples
#' NEED TO INSERT
#' 
#' @export

calculateSaturationVapourPressure <- function(temp){
      
      satVapPres <- .611 * exp((17.27 * temp) / (temp + 237.3))
      
      return(satVapPres)
      
}

#' @title calculateSaturationVapourPressure.
#'
#' @description
#' \code{calculateSaturationVapourPressure} function to return the specific heat of moist air
#' 
#' @details
#' function to return the specific heat of moist air
#' variable in FAO 56 3-10
#'
#' @references see aWhere API Documentation (insert HTTP address)
#' 
#' @return double
#'
#' @examples
#' NEED TO INSERT
#' 
#' @export

calculateSpecificHeatAir <- function() {
      
      return(1.013)
}

#' @title calculateAerodynamicResistance.
#'
#' @description
#' \code{calculateAerodynamicResistance} function to calculate aerodynamic resistance for a given crop
#' 
#' @details
#' function to calculate aerodynamic resistance for a given crop
#' takes crop height (m) and wind speed as inputs
#' EQN 4 in FAO 56
#' 
#' zW height of wind measurements [m],
#' zH height of humidity measurements [m],
#' zero plane displacement height [m],
#' zOM roughness length governing momentum transfer [m],
#' zOH roughness length governing transfer of heat and vapour [m],
#' k von Karman's constant, 0.41 [-],
#' z wind speed at height z [m s-1].
#'
#' @references see aWhere API Documentation (insert HTTP address)
#'
#' @param - windspeed: In m/s measured at 2m - double
#' @param - cropHeight: in m
#' 
#' @return double
#'
#' @examples
#' NEED TO INSERT
#' 
#' @export

calculateAerodynamicResistance <- function(windSpeed, cropHeight) {
      
      zeroPlaneDisHeight <- calculateZeroPlaneDisplacement(cropHeight)
      zOM <- calculateRoughnessLength(cropHeight)
      zOH <- .1 * zOM
      
      #these are hardcoded and are due to weather weather instrumentation is
      #located relative to ground
      zW <- 10
      zH <- 2
      vonKarman <- .41
      
      aeroResis <- (log((zW - zeroPlaneDisHeight)/zOM) * 
                          log((zH - zeroPlaneDisHeight)/zOH)) / (vonKarman^2 * windSpeed)
      
      return(aeroResis)
}

#' @title calculateZeroPlaneDisplacement.
#'
#' @description
#' \code{calculateZeroPlaneDisplacement} function to calculate the zero plane displacement  
#' 
#' @details
#' function to calculate the zero plane displacement  for aerodynamic resistance
#' calculations for a given crop takes crop height (m).  Values from FAO 56
#'
#' @references see aWhere API Documentation (insert HTTP address)
#'
#' @param - cropHeight: In m - double
#' 
#' @return double
#'
#' @examples
#' NEED TO INSERT
#' 
#' @export

calculateZeroPlaneDisplacement <- function(cropHeight) {
      
      zeroPlaneDisHeight <- (2/3) * cropHeight
      
      return(zeroPlaneDisHeight)
}

#' @title calculateRoughnessLength.
#'
#' @description
#' \code{calculateRoughnessLength} function to calculate the roughness length 
#' 
#' @details
#' function to calculate the roughness length for aerodynamic resistance calculations for a given crop
#' takes crop height (m).  Note that I am using the values from FAO 56 instead of the CART model paper
#' EQN 4 in FAO 56
#'
#' @references see aWhere API Documentation (insert HTTP address)
#'
#' @param - cropHeight: In m - double
#' 
#' @return double
#'
#' @examples
#' NEED TO INSERT
#' 
#' @export

calculateRoughnessLength <- function(cropHeight) {
      
      roughnessLength <- .123 * cropHeight
      
      return(roughnessLength)
}

#' @title calculateActualVapourPressure.
#'
#' @description
#' \code{calculateActualVapourPressure} function to calculate the actual vapour pressure
#' 
#' @details
#' function to calculate the actual vapour pressure from min/max RH and temp
#' eqn 15 in FAO 56
#'
#' @references see aWhere API Documentation (insert HTTP address)
#'
#' @param - rhMin: double
#' @param - rhMax: double
#' @param - tMin: double
#' @param - tMax: double
#' 
#' @return double
#'
#' @examples
#' NEED TO INSERT
#' 
#' @export

calculateActualVapourPressure <- function(rhMin,rhMax,tMin,tMax) {
      
      satVapPresL <- calculateSaturationVapourPressure(tMin)
      satVapPresH <- calculateSaturationVapourPressure(tMax)
      
      actVapPres <- ((satVapPresL *(rhMax/100)) + (satVapPresH * (rhMin/100)))/2
 
      return(actVapPres)     
}

#' @title calculatePsychrometricConstant.
#'
#' @description
#' \code{calculatePsychrometricConstant} function to calculate the psychometric contstant
#' 
#' @details
#' function to calculate the physchometric constant based on on atmospheric
#' pressure/elevation
#' eqn 3-10 in FAO 56
#'
#' @references see aWhere API Documentation (insert HTTP address)
#'
#' @param - elevation: in m - double
#' @param - elevAtRefLevel: in m - double
#' @param - temp: in C double
#' 
#' @return double
#'
#' @examples
#' NEED TO INSERT
#' 
#' @export

calculatePsychrometricConstant <- function(elevation,elevAtRefLevel,temp){
      
      ratioMolWeightWetDryAir <- .622
      specHeat <- calculateSpecificHeatAir()
      atmPres <- calculateAtmosphericPressure(elevation,elevAtRefLevel,temp)
      latentHeatVap <- calculateLatentHeatVap() 
      
      psychoConst <- ((specHeat * atmPres)/(ratioMolWeightWetDryAir * latentHeatVap)) * 10^(-3)
      
      return(psychoConst)
}

#' @title calculateDewPointDepression.
#'
#' @description
#' \code{calculateDewPointDepression} function to calculate the dew point depression
#' 
#' @details
#' function to calculate Dew Point Depression
#' definition from Kim, 2002 paper
#'
#' @references see aWhere API Documentation (insert HTTP address)
#'
#' @param - airTemp: in C - double
#' @param - RH: 1-100 - double
#' @param - temp: in C double
#' 
#' @return double
#'
#' @examples
#' NEED TO INSERT
#' 
#' @export

calculateDewPointDepression <- function(airTemp,RH,temp){
      
      actualVapPres <- calculateActualVapourPressure(RH,RH,temp,temp)
      dewPointTemp <- calculateDewPointTemp(actualVapPres)
      dpd <- airTemp - dewPointTemp
 
      return(dpd)     
}

#' @title calculateRosenbergWind.
#'
#' @description
#' \code{calculateRosenbergWind} function to calculate adjusted windspeed
#' 
#' @details
#' function to calculate Rosenberg's Adjusted windspeed for the CART model info
#' comes http://agsys.cra-cin.it/tools/leafwetness/help/ 
#' paper is Model to Enhance Site-Specific Estimation of Leaf Wetness Duration by Kim et al, 2002
#'
#' @references see aWhere API Documentation (insert HTTP address)
#'
#' @param - wind: in m/s - double
#' @param - cropHeight: in m - double
#' 
#' @return double
#'
#' @examples
#' NEED TO INSERT
#' 
#' @export

calculateRosenbergWind <- function(wind, cropHeight) {
      
      zpH <- calculateZeroPlaneDisplacement(cropHeight)
      RL  <- calculateRoughnessLength(cropHeight)
      
      RWS <- wind * ((log(cropHeight - zpH) - log(RL))/(log(10-zpH) - log(RL)))
      
      return(RWS)
}

#' @title calculateLandsbergWind.
#'
#' @description
#' \code{calculateLandsbergWind} function to calculate adjusted windspeed
#' 
#' @details
#' function to calculate Lansberg's Adjusted windspeed
#'
#' @references see aWhere API Documentation (insert HTTP address)
#'
#' @param - wind: in m/s - double
#' @param - cropHeight: in m - double
#' @param - canopyHeight: in m - double
#' 
#' @return double
#'
#' @examples
#' NEED TO INSERT
#' 
#' @export

calculateLandsbergWind <- function(wind, cropHeight, canopyHeight) {
      
      LWS <- wind * (1 + (1/3) * ((1-canopyHeight)/cropHeight))^-2
      
      return(LWS)
}

#' @title calculateBoundaryLayerResistance.
#'
#' @description
#' \code{calculateBoundaryLayerResistance} function to calculate boundary layer resistance
#' 
#' @details
#' function to calculate Boundary Layer Resistance for LWR model info
#' comes http://agsys.cra-cin.it/tools/leafwetness/help/ 
#'
#' @references see aWhere API Documentation (insert HTTP address)
#'
#' @param - dimMockLeaf: double
#' @param - wind: in m/s - double
#' 
#' @return double
#'
#' @examples
#' NEED TO INSERT
#' 
#' @export

calculateBoundaryLayerResistance <- function(dimMockLeaf, wind) {
      
      rBL <- 307*sqrt(dimMockLeaf/wind)
      
      return(rBL)
}