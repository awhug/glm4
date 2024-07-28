### glm4 methods ----

#' @inherit stats::family
#' @export
#' @keywords internal
family.glm4 <- function (object, ...) object$family

#' @inherit stats::logLik
#' @export
#' @keywords internal
logLik.glm4 <- function (object, ...){
	if (!missing(...))
		warning("extra arguments discarded")
	fam <- family.glm4(object)$family
	p <- object$rank[1]
	if (fam %in% c("gaussian", "Gamma", "inverse.gaussian"))
		p <- p + 1
	val <- p - object$aic/2
	attr(val, "nobs") <- sum(!is.na(object$residuals))
	attr(val, "df") <- p
	class(val) <- "logLik"
	val
}

#' @inherit stats::predict.glm
#' @export
#' @keywords internal
predict.glm4 <- function(...) stats:::predict.glm(...)

#' @inherit stats::residuals.glm
#' @export
#' @keywords internal
residuals.glm4 <- function(object, ...) MatrixModels:::residuals(object = object$glm4_fit, ...)

#' @inherit stats::residuals.glm
#' @export
#' @keywords internal
resid.glm4 <- residuals.glm4

#' @inherit stats::vcov
#' @export
#' @keywords internal
vcov.glm4 <- function (object, complete = TRUE, ...) stats:::vcov.summary.glm(summary.glm4(object, ...), complete = complete)

#' @inherit stats::confint
#' @export
#' @keywords internal
confint.glm4 <- function (object, parm, level = 0.95, ...)
{
	cf <- coef(object)
	pnames <- names(cf)
	if (missing(parm))
		parm <- pnames
	else if (is.numeric(parm))
		parm <- pnames[parm]
	a <- (1 - level)/2
	a <- c(a, 1 - a)
	pct <- stats:::format.perc(a, 3)
	fac <- qnorm(a)
	ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, pct))
	ses <- sqrt(Matrix::diag(vcov(object)))[parm]
	ci[] <- cf[parm] + ses %o% fac
	ci
}

#' @inherit stats::df.residual
#' @export
#' @keywords internal
df.residual.glm4 <- function(object){
	if (!is(object, "glpModel")){
		object <- object$glm4_fit
	}
	fit <- object
	rank <- as.numeric(rank.glm4(object))
	Xdims <- dim(object@pred@X)
	nobs <- Xdims[1]
	n.ok <- nobs - sum(object@resp@weights==0)
	resdf <- n.ok - rank
	return(resdf)
}

# Convenience function to get the rank

#' @export
#' @keywords internal
rank.glm4 <- function(object) Matrix::rankMatrix(object@pred@X, method = "qr.R")

#' @export
#' @keywords internal
print.glm4 <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
	cat("\nCall:  ",
			paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
	if(length(coef(x))) {
		cat("Coefficients")
		if(is.character(co <- x$contrasts))
			cat("  [contrasts: ", apply(cbind(names(co),co), 1L, paste, collapse = "="), "]")
		cat(":\n")
		print.default(
			x = format(x$coefficients, digits = digits),
			print.gap = 2,
			quote = FALSE)
	} else cat("No coefficients\n\n")
	cat("\nDegrees of Freedom:", x$df.null, "Total (i.e. Null); ",
			x$df.residual, "Residual\n")
	if(nzchar(mess <- naprint(x$na.action))) cat("  (",mess, ")\n", sep = "")
	cat("Null Deviance:	   ",	format(signif(x$null.deviance, digits)),
			"\nResidual Deviance:", format(signif(x$deviance, digits)),
			"\tAIC:", format(signif(x$aic, digits)))
	cat("\n")
	invisible(x)
}

#' @export
#' @keywords internal
weights.glm4 <- function(object, type = c("prior", "working"), ...){
	type <- match.arg(type)
	res <- if(type == "prior") object$prior.weights else object$weights
	if(is.null(object$na.action)) res
	else naresid(object$na.action, res)
}

#' @export
#' @keywords internal
model.matrix.glm4 <- function(object, Matrix = FALSE){
	if (!is(object, "glm4")){
		stop("Model has not been fit using glm4/MatrixModels")
	}
	mm <- object$model
	if (!Matrix){
		assign_vals <- attr(mm, "assign")
		mm <- as.matrix(mm)
		attr(mm, "assign") <- assign_vals
	}
	return(mm)
}
