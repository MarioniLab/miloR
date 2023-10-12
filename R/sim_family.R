#' sim_family
#'
#' Simulated counts data from a series of simulated family trees
#'
#' Data are simulated counts from 30 families and includes X and Z design matrices,
#' as well as a single large kinship matrix. Kinships between family members are
#' dictated by the simulated family, i.e. sibs=0.5, parent-sib=0.5, sib-grandparent=0.25, etc.
#' These kinships, along with 2 other random effects, are used to induce a defined covariance
#' between simulated obserations as such:
#'
#' Z:= random effect design matrix, n X q
#' G:= matrix of variance components, including kinship matrix
#'
#' LL^T = Chol(ZGZ^T) := the Cholesky decomposition of the random effect contribution
#' to the sample covariance
#' Ysim:= simulated means based on exp(offset + Xbeta + Zb)
#' Y = LYsim := simulated means with defined covariance
#'
#'
#' @docType data
#' @usage data(sim_family)
#'
#' @format A list containing a \code{data.frame} in the "DF" slot containing the
#' mean counts and meta-data, and a \code{matrix} containing the kinship matrix
#' across all families in the "IBD" slot.
#'
#' @keywords datasets
#'
#'
#' @examples
#' NULL
#'
#' @name sim_family
#' @aliases sim_family
NULL
