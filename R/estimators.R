# ============================================================
# Four exported causal estimators.  Each takes a tidy data
# frame / vectors and returns a named list(beta, se, pvalue).
# ============================================================


#helpers (internal)
.covar_str <- function(covar_names) {
  if (length(covar_names) == 0) return("")
  paste0(" + ", paste(covar_names, collapse = " + "))
}

.extract_coef <- function(fit, term) {
  sm <- summary(fit)$coefficients
  list(b  = coef(fit)[term],
       se = sm[term, 2],
       p  = sm[term, 4])
}


#1. DIRECT

#' DIRECT estimator: OLS with instrument and negative-control as covariates
#'
#' Regresses Y on Z plus the genetic instrument G, the negative-control W,
#' and any additional covariates.  This is a "naïve" adjustment that uses
#' whatever observables are available but does NOT correct for unmeasured
#' confounding via a ratio or IV approach.
#'
#' @param y       Numeric outcome vector (length n).
#' @param Z       Numeric exposure vector (length n), assumed pre-scaled.
#' @param g       Numeric instrument vector (length n).
#' @param w       Numeric negative-control vector (length n).
#' @param covars  Optional data frame of additional covariates (n rows).
#'
#' @return Named list: \code{beta}, \code{se}, \code{pvalue}.
#' @export
#'
#' @examples
#' set.seed(1)
#' n   <- 200
#' dat <- generate_toy_data(n = n, seed = 1)  # internal helper
#' fit_direct(dat$Y[, 1], dat$Z, dat$G[, 1], dat$W[, 1], dat$synthetic_data)
fit_direct <- function(y, Z, g, w, covars = NULL) {
  NA_result <- list(beta = NA_real_, se = NA_real_, pvalue = NA_real_)
  cnames <- if (!is.null(covars)) names(covars) else character(0)
  d      <- cbind(data.frame(y = y, Z = Z, g = g, w = w), covars)
  fml    <- as.formula(paste0("y ~ Z + g + w", .covar_str(cnames)))

  fit <- tryCatch(lm(fml, data = d), error = function(e) NULL)
  if (is.null(fit)) return(NA_result)

  sm <- summary(fit)$coefficients
  if (!"Z" %in% rownames(sm)) return(NA_result)

  list(
    beta   = as.numeric(coef(fit)["Z"]),
    se     = as.numeric(sm["Z", 2]),
    pvalue = as.numeric(sm["Z", 4])
  )
}

#2. COCA

#' COCA estimator: Negative-Control Outcome Correction via ratio
#'
#' Implements the Correlated Outcome Control Approach (COCA).  Fits
#' \code{w ~ y + Z + covars} and recovers the causal effect as
#' \eqn{\hat\beta = -\hat\beta_Z / \hat\beta_Y}.  Standard errors are
#' obtained via the delta method.
#'
#' The negative-control W should be an outcome that shares the same
#' unmeasured confounders as Y but has no direct causal path from Z.
#'
#' @param y       Numeric primary outcome vector (length n).
#' @param Z       Numeric exposure vector (length n).
#' @param w       Numeric negative-control outcome vector (length n).
#'                Recommended: pass \code{rowMeans(W_matrix)} for stability.
#' @param covars  Optional data frame of additional covariates (n rows).
#' @param ratio_cap  Maximum absolute value of the ratio estimate before the
#'                   result is flagged as unstable and \code{NA} is returned.
#'                   Default 10.
#' @param se_cap     Maximum SE before flagging as unstable. Default 5.
#'
#' @return Named list: \code{beta}, \code{se}, \code{pvalue}.
#'         Returns \code{list(beta=NA, se=NA, pvalue=NA)} if estimation is
#'         unstable (near-zero \eqn{\hat\beta_Y} or extreme ratio).
#' @export

fit_coca <- function(y, Z, w, covars = NULL,
                     ratio_cap = 10, se_cap = 5) {
  NA_result <- list(beta = NA_real_, se = NA_real_, pvalue = NA_real_)
  cnames    <- if (!is.null(covars)) names(covars) else character(0)
  d         <- cbind(data.frame(w = w, y = y, Z = Z), covars)
  fml       <- as.formula(paste0("w ~ y + Z", .covar_str(cnames)))

  fit <- tryCatch(lm(fml, data = d), error = function(e) NULL)
  if (is.null(fit)) return(NA_result)

  b  <- coef(fit)
  bZ <- b["Z"]
  bY <- b["y"]

  if (is.na(bY) || abs(bY) < 1e-8) return(NA_result)

  bhat <- -bZ / bY
  V    <- vcov(fit)[c("Z", "y"), c("Z", "y")]
  grad <- c(-1 / bY, bZ / bY^2)
  se_h <- sqrt(as.numeric(t(grad) %*% V %*% grad))

  if (abs(bhat) > ratio_cap || is.nan(se_h) || se_h > se_cap) return(NA_result)

  list(
    beta   = bhat,
    se     = se_h,
    pvalue = 2 * (1 - pnorm(abs(bhat / se_h)))
  )
}


#3. IV2SLS

#' IV2SLS estimator: Two-Stage Least Squares with genetic instrument
#'
#' Uses the genetic instrument G to instrument for the exposure Z,
#' controlling for the negative-control W and any additional covariates.
#' Requires \pkg{AER}.
#'
#' A weak-instrument check (first-stage F < 10) is applied; if the
#' instrument is weak the function returns \code{NA}.
#'
#' @param y       Numeric outcome vector (length n).
#' @param Z       Numeric exposure vector (length n).
#' @param g       Numeric instrument vector (length n).
#' @param w       Numeric negative-control vector (length n).
#' @param covars  Optional data frame of additional covariates (n rows).
#' @param min_f   Minimum acceptable first-stage F-statistic. Default 10.
#'
#' @return Named list: \code{beta}, \code{se}, \code{pvalue}.
#' @export
fit_iv2sls <- function(y, Z, g, w, covars = NULL, min_f = 10) {
  NA_result <- list(beta = NA_real_, se = NA_real_, pvalue = NA_real_)
  cnames    <- if (!is.null(covars)) names(covars) else character(0)
  cs        <- .covar_str(cnames)

  # First-stage weak-instrument check
  d_fs  <- cbind(data.frame(Z = Z, g = g, w = w), covars)
  fml_fs <- as.formula(paste0("Z ~ g + w", cs))
  fs    <- tryCatch(lm(fml_fs, data = d_fs), error = function(e) NULL)
  if (is.null(fs)) return(NA_result)
  Fst   <- summary(fs)$fstatistic[1]
  if (is.na(Fst) || Fst < min_f) return(NA_result)

  # 2SLS
  d_iv  <- cbind(data.frame(y = y, Z = Z, G_inst = g, w = w), covars)
  fml_2sls <- as.formula(
    paste0("y ~ Z + w", cs, " | G_inst + w", cs)
  )
  fit <- tryCatch(
    AER::ivreg(fml_2sls, data = d_iv),
    error = function(e) NULL
  )
  if (!"Z" %in% names(coef(fit))) return(NA_result)
  list(
    beta   = as.numeric(coef(fit)["Z"]),
    se     = as.numeric(sm["Z", 2]),
    pvalue = as.numeric(sm["Z", 4])
  )
}


#4. PGC

#' PGC estimator: Proxy G-Component Correction
#'
#' A three-step estimator that:
#' \enumerate{
#'   \item Removes the genetic (G) component from Z to isolate the
#'         U-driven residual.
#'   \item Bridges the negative-control W onto that U-residual to
#'         construct \eqn{\hat W}, a proxy for unmeasured confounding.
#'   \item Includes \eqn{\hat W} in the final outcome regression to
#'         partially absorb confounding bias.
#' }
#'
#' @param y       Numeric outcome vector (length n).
#' @param Z       Numeric exposure vector (length n).
#' @param g       Numeric instrument vector (length n).
#' @param w       Numeric negative-control vector (length n).
#'                Pass \code{rowMeans(W_matrix)} for stability.
#' @param covars  Optional data frame of additional covariates (n rows).
#'
#' @return Named list: \code{beta}, \code{se}, \code{pvalue}.
#' @export
fit_pgc <- function(y, Z, g, w, covars = NULL) {
  NA_result <- list(beta = NA_real_, se = NA_real_, pvalue = NA_real_)
  cnames    <- if (!is.null(covars)) names(covars) else character(0)
  cs        <- .covar_str(cnames)

  # Step 1: residualise Z on G -> U-driven residual
  d_r       <- cbind(data.frame(Zc = Z, g = g), covars)
  fit_resid <- tryCatch(
    lm(as.formula(paste0("Zc ~ g", cs)), data = d_r),
    error = function(e) NULL
  )
  if (is.null(fit_resid)) return(NA_result)
  Z_resid <- residuals(fit_resid)

  # Step 2: bridge W on the U-component
  d_b   <- cbind(data.frame(w = w, Z_resid = Z_resid), covars)
  fit_b <- tryCatch(
    lm(as.formula(paste0("w ~ Z_resid", cs)), data = d_b),
    error = function(e) NULL
  )
  if (is.null(fit_b)) return(NA_result)
  W_hat <- fitted(fit_b)

  # Step 3: final outcome regression
  d_f   <- cbind(data.frame(y = y, Z = Z, W_hat = W_hat), covars)
  fit_f <- tryCatch(
    lm(as.formula(paste0("y ~ Z + W_hat", cs)), data = d_f),
    error = function(e) NULL
  )
  if (!"Z" %in% names(coef(fit))) return(NA_result)
  list(
    beta   = as.numeric(coef(fit)["Z"]),
    se     = as.numeric(sm["Z", 2]),
    pvalue = as.numeric(sm["Z", 4])
  )
}
