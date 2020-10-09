#' @title sim_discrete
#' Simulated discrete groups data
#'
#' Data are simulated single-cells in 4 distinct groups of cells. Cells in each
#' group are assigned to 1 of 2 conditions: \emph{A} or \emph{B}. Specifically,
#' the cells in block 1 are highly abundant in the \emph{A} condition, whilst
#' cells in block 4 are most abundant in condition \emph{B}.
#'
#' @docType data
#'
#' @usage data(sim_discrete)
#'
#' @format A list containing a \code{\linkS4class{Milo}} object in the "mylo" slot,
#' and a \code{data.frame} containing experimental meta-data in the "meta" slot.
#'
#' @keywords datasets
#'
#'
#' @examples
#' load(sim_discrete)
#' print(sim_discrete$mylo)
#'
#' head(sim_discrete$meta)
NULL
