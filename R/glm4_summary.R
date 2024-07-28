#' Summarising Generalized Linear Models Using Sparse Matrices
#'
#' @description
#' Generates a summary of the `glm4` object to evaluate coefficients, standard errors, model fit criteria, etc.
#'
#' @usage
#' ## S3 method for class 'glm4'
#' summary(object, p.adjust = NULL, dispersion = NULL, correlation = FALSE, symbolic.cor = FALSE, ...)
#'
#' @param object an object of class "glm4"
#' @param p.adjust returns p-values adjusted using one of several methods implemented in `stats::p.adjust`. Defaults to `NULL` for no adjustment, consistent with `stats::glm`.
#' @param dispersion the dispersion parameter for the family used. Either a single numerical value or NULL (the default), when it is inferred from object (see `stats::summary.glm()` details).
#' @param correlation logical; if `TRUE`, the correlation matrix of the estimated parameters is returned and printed.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' This function is designed to resemble `stats::summary.glm()` as closely as possible. It calculates the (un)scaled variance-covariance matrix from a `MatrixModels::glm4()` object using `Matrix::chol2inv()` and produces a coefficient table for easy inspection of model parameters.
#'
#' @returns A list object of class `c("summary.glm", "summary.glm4")`. See `stats::summary.glm()` for more details on returned components.
#' @export summary.glm4
#' @export

summary.glm4 <- function(
		object, p.adjust = NULL, dispersion = NULL, correlation = FALSE, symbolic.cor = FALSE,
		...)
{
	fit <- object$glm4_fit
	est.disp <- FALSE
	df.r <- df.residual.glm4(fit)
	mod.residuals <- MatrixModels::residuals(fit, type = "working")
	if (is.null(dispersion))
		dispersion <- if (object$family$family %in% c("poisson",
																									"binomial"))
			1
	else if (df.r > 0) {
		est.disp <- TRUE
		if (any(object$weights == 0))
			warning("observations with zero weight not used for calculating dispersion")
		sum((object$weights * mod.residuals^2)[object$weights >
																					 	0])/df.r
	}
	else {
		est.disp <- TRUE
		NaN
	}

	coef.p <- MatrixModels::coef(fit)
	aliased <- is.na(coef.p)
	p <- rank.glm4(fit)
	if (p > 0){
		p1 <- 1L:p
		fac <- fit@pred@fac
		# When sparse = FALSE, fac returns an unnamed dpoMatrix.
		#	When sparse = TRUE, it returns a named dsCMatrix.
		# Need to coerce and name dims
		if (is(fac, "Cholesky")){
			dnms <- fac@Dimnames
			fac <- as(fac, "dtrMatrix")
			rownames(fac) <- colnames(fac) <- dnms[[1]]
		}
		covmat.unscaled <- Matrix::chol2inv(fac)
		covmat <- dispersion * covmat.unscaled
		var.cf <- Matrix::diag(covmat)
		s.err <- sqrt(var.cf)
		tvalue <- coef.p/s.err
		dn <- c("Estimate", "Std. Error")
		if (!est.disp) {
			pvalue <- 2 * pnorm(-abs(tvalue))
			coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
			dimnames(coef.table) <- list(names(coef.p), c(dn,
																										"z value", "Pr(>|z|)"))
		}
		else if (df.r > 0) {
			pvalue <- 2 * pt(-abs(tvalue), df.r)
			coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
			dimnames(coef.table) <- list(names(coef.p), c(dn,
																										"t value", "Pr(>|t|)"))
		}
		else {
			coef.table <- cbind(coef.p, NaN, NaN, NaN)
			dimnames(coef.table) <- list(names(coef.p), c(dn,
																										"t value", "Pr(>|t|)"))
		}
		# Adjust the p.value
		if (!is.null(p.adjust)){
			coef.table[,4] <- stats::p.adjust(coef.table[,4], method = p.adjust)
		}
		df.f <- ncol(fac)
	}
	else {
		coef.table <- matrix(nrow=0L, ncol=4L)
		dimnames(coef.table) <- list(NULL, c("Estimate", "Std. Error",
																				 "t value", "Pr(>|t|)"))
		covmat.unscaled <- covmat <- matrix(nrow=0L, ncol=0L)
		df.f <- length(aliased)
	}

	## return answer

	## these need not all exist, e.g. na.action.
	keep <- match(c("call","terms","family","deviance", "aic",
									"contrasts", "df.residual","null.deviance","df.null",
									"iter", "na.action"), names(object), 0L)
	ans <- c(object[keep],
					 list(deviance.resid = MatrixModels::residuals(fit, type = "deviance"),
					 		 coefficients = coef.table,
					 		 aliased = aliased,
					 		 dispersion = dispersion,
					 		 df = c(object$rank, df.r, df.f),
					 		 cov.unscaled = covmat.unscaled,
					 		 cov.scaled = covmat))

	if(correlation && p > 0) {
		dd <- sqrt(diag(covmat.unscaled))
		ans$correlation <-
			covmat.unscaled/outer(dd,dd)
		ans$symbolic.cor <- symbolic.cor
	}

	# Setting the class as summary.glm allows generic methods
	# including (importantly) print.summary.glm
	class(ans) <- c("summary.glm", "summary.glm4")
	return(ans)
}
