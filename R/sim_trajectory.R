#' @title sim_trajectory
#' Simulated linear trajectory data
#'
#' Data are simulated single-cells along a single linear trajectory. Cells are
#' simulated from 5 groups, and assigned to 1 of 2 conditions; \emph{A} or \empj{B}.
#' Data were generated using in the \code{\link[dyntoy]{simulate_linear_trajectory}}
#' function in the \code{dyntoy} package.
#'
#' @docType data
#' @usage data(sim_trajectory)
#'
#' @format A list containing a \code{\linkS4class{Milo}} object in the "mylo" slot,
#' and a \code{data.frame} containing experimental meta-data in the "meta" slot.
#'
#' @keywords datasets
#'
#' @references https://github.com/dynverse/dyntoy
#'
#' @examples
#' data(sim_trajectory)
#' print(sim_trajectory$mylo)
#'
#' head(sim_trajectory$meta)
NULL
