#' sim_nbglmm
#'
#' Simulated counts data from a NB-GLMM for a single trait
#'
#' Data are simulated counts from 50 samples in a single data frame, from which the
#' X and Z design matrices, can be constructed (see examples). There are 2 random effects and 2 fixed
#' effect variables used to simulate the count trait.
#'
#'
#' @docType data
#' @usage data(sim_nbglmm)
#'
#' @format A \code{data.frame} \emph{sim_nbglmm} containing the following columns:
#' \describe{
#' \item{\code{Mean:}}{\code{numeric} containing the base mean computed as the linear combination of the
#' simulated fixed and random effect weights multiplied by their respective weight matrices.}
#' \item{\code{Mean.Count:}}{\code{numeric} containing the integer count values randomly sampled from a negative
#' binomail distribution with mean = \emph{Mean} and dispersion = \emph{r}}
#' \item{\code{r:}}{\code{numeric} containing the dispersion value used to simulate the integer counts in
#' \emph{Mean.Count}.}
#' \item{\code{Intercept:}}{\code{numeric} of all 1s which can be used to set the intercept term in the X design
#' matrix.}
#' \item{\code{FE1:}}{\code{numeric} a binary fixed effect variable taking on values [0, 1]}
#' \item{\code{FE2:}}{\code{numeric} a continuous fixed effect variables}
#' \item{\code{RE1:}}{\code{numeric} a random effect variable with 10 levels}
#' \item{\code{RE2:}}{\code{numeric} a random effect variable with 7 levels}
#' }
#'
#' @keywords datasets
#'
#' @examples
#' data(sim_nbglmm)
#' head(sim_nbglmm)
#'
#' @name sim_nbglmm
#' @aliases sim_nbglmm
NULL
