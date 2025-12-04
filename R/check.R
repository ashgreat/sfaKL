#' Lightweight diagnostic check for an sfaKL fit
#'
#' Validates positive definiteness of key matrices and reports condition numbers.
#'
#' @param object An object of class sfaKL
#' @return A list with PD flags and condition numbers
#' @export
sfaKL_check <- function(object) {
    stopifnot(inherits(object, "sfaKL"))

    mats <- try(get_matrices(object$raw_par, J = object$J, M = object$M), silent = TRUE)
    if (inherits(mats, "try-error")) {
        return(list(ok = FALSE, error = "Failed to build matrices from parameters"))
    }

    is_pd <- function(M) {
        eigenvalues <- try(eigen(M, symmetric = TRUE, only.values = TRUE)$values, silent = TRUE)
        if (inherits(eigenvalues, "try-error")) return(FALSE)
        all(eigenvalues > 0)
    }

    list(
        ok = TRUE,
        Theta_pd = is_pd(mats$Theta),
        Delta_pd = is_pd(mats$Delta),
        Sigma_pd = is_pd(mats$Sigma),
        Omega_pd = is_pd(mats$Omega),
        kappa_Theta = kappa(mats$Theta),
        kappa_Delta = kappa(mats$Delta)
    )
}
