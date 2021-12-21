#' COVID-19 case numbers and environmental factors.
#'
#' A dataset containing the COVID-19 case numbers and 
#' 7 environmental factors, for 7 counties in New York
#' State between March 1, 2020 and September 30, 2021.
#' 
#' @format A data frame with 3600 rows and 32 variables:
#' \describe{
#'   \item{time}{date of the observation}
#'   \item{county}{county of the observation}
#'   \item{temp}{temperature, in Fahrenheit}
#'   \item{dew}{dew point, in Fahrenheit}
#'   \item{wind}{wind speed, in miles per hour}
#'   \item{precipitation}{precipitation, in inches}
#'   \item{humidity}{humidity, in percentage}
#'   \item{pm25}{fine particles with an aerodynamic 
#'   diameter of 2.5microgram or less, 
#'   in microgram per cubic meter}
#'   \item{ozone}{ozone, in microgram per cubic meter}
#'   \item{case}{COVID-19 case number}
#'   \item{case_avg}{COVID-19 case number weekly average, 
#'   between each day and following 6 days}
#'   \item{case_avg_next_n}{COVID-19 case number weekly average,
#'   for n day after current date}
#' }
#' @source \url{http://www.nrcc.cornell.edu/}
#' @source \url{https://www.epa.gov/}
#' @source \url{https://data.ny.gov/}
"covid.environment"
