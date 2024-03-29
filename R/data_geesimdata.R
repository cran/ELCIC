#' Data simulated for model selection under GEE framework without missingness
#'
#'
#' @docType data
#'
#' @usage data(geesimdata)
#'
#' @format An object of class \code{"list"}
#' \describe{
#'  \item{y}{The outcomes generated from Poisson distribution with three repeated measurements from each subject}
#'  \item{x}{A covariate matrix, of which the first column are all ones and rest columns contain normally distributed. Two are time-dependent variables, and one is time-independent variable.}
#' }
#' @references This data set was artificially created for the ELCIC package.
#' @keywords datasets
#' @examples
#'
#' data(geesimdata)
#' geesimdata$y
#'
"geesimdata"
