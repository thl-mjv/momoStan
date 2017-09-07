#' Finnish mortality data
#'
#' A dataset containing Finnish mortality data in five agegroups 
#'
#' @format A data frame with ten variables and 2040 (5*340) rows
#' \describe{
#'   \item{date}{date of the Monday of the week}
#'   \item{n}{number of deaths}
#'   \item{pop}{approximate population}
#'   \item{TEMP}{weekly mean temperature}
#'   \item{TMIN}{weekly mean daily max temperature}
#'   \item{TMAX}{weekly mean daily min temperature}
#'   \item{DEWP}{weekly mean dew point}
#'   \item{SLP}{weekly mean sea level air pressure}
#'   \item{InfA}{number of Influenza notifications to the National Infectious Disease register}
#'   \item{age}{age group}
#' }
#' @source VRK, TTR, NOAA
"fimomodata"
