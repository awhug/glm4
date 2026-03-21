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
	#pct <- stats:::format.perc(a, 3)
	pct <- paste(format(100 * a, trim = TRUE, scientific = FALSE, digits = 3),"%")
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

#' @inherit stats::case.names
#' @export
#' @keywords internal
case.names.glm4 <- function(object, full = FALSE, ...) {
	w <- weights(object)
	dn <- names(object$residuals)
	if (full || is.null(w)) dn else dn[w != 0]
}

#' @inherit stats::variable.names
#' @export
#' @keywords internal
variable.names.glm4 <- function(object, full = FALSE, ...) {
	if (full) colnames(model.matrix(object)) else names(coef(object))
}

#' @keywords internal
dispersion.glm4 <- function(object, dispersion = NULL) {
	if (!is.null(dispersion)) return(dispersion)
	if (object$family$family %in% c("poisson", "binomial")) return(1)
	df.r <- object$df.residual
	if (df.r > 0) {
		mod.residuals <- MatrixModels:::residuals(object$glm4_fit, type = "working")
		if (any(object$weights == 0))
			warning("observations with zero weight not used for calculating dispersion")
		sum((object$weights * mod.residuals^2)[object$weights > 0]) / df.r
	} else {
		NaN
	}
}

#' @inherit stats::nobs
#' @export
#' @keywords internal
nobs.glm4 <- function(object, ...) sum(!is.na(object$residuals))

#' @inherit stats::deviance
#' @export
#' @keywords internal
deviance.glm4 <- function(object, ...) object$deviance

#' @inherit stats::formula
#' @export
#' @keywords internal
formula.glm4 <- function(x, ...) x$formula

#' @inherit stats::sigma
#' @export
#' @keywords internal
sigma.glm4 <- function(object, ...) sqrt(dispersion.glm4(object))

#' @inherit stats::extractAIC
#' @export
#' @keywords internal
extractAIC.glm4 <- function(fit, scale = 0, k = 2, ...) {
	edf <- fit$rank
	c(edf, fit$aic + (k - 2) * edf)
}

#' @inherit stats::model.frame
#' @export
#' @keywords internal
model.frame.glm4 <- function(formula, ...) {
	model.frame(formula$terms, formula$data)
}

#' @export
#' @keywords internal
labels.glm4 <- function(object, ...) attr(object$terms, "term.labels")

### influence methods ----

# All three measures share the hat values h_i = w_i * x_i^T (X^TWX)^{-1} x_i,
# computed via the Cholesky factor in glm4_fit

#' @inherit stats::hatvalues
#' @param batch_size integer; number of rows processed per batch when the model matrix is sparse. Reduce if memory is limited.
#' @param verbose logical; if `TRUE`, a progress bar is printed during the batched sparse computation.
#' @export
#' @keywords internal
hatvalues.glm4 <- function(model, batch_size = 1000L, verbose = FALSE, ...) {
	fac <- model$glm4_fit@pred@fac
	X <- model$glm4_fit@pred@X
	w <- model$weights
	n <- length(w)

	if (is(fac, "CHMfactor")) {
		# Batched to avoid a p*n dense intermediate (~GB for large models).
		h <- numeric(n)
		batches <- seq(1L, n, by = batch_size)
		if (verbose) {
			pb <- utils::txtProgressBar(min = 0L, max = length(batches), style = 3)
			on.exit(close(pb))
		}
		# Each batch solves (X^TWX)^{-1} X_batch^T (p×B), then extracts the diagonal
		# of X_batch %*% V_batch via rowSums(X_batch * t(V_batch)).
		for (i in seq_along(batches)) {
			s <- batches[i]
			idx <- s:min(s + batch_size - 1L, n)
			X_b <- X[idx, , drop = FALSE]
			V_b <- Matrix::solve(fac, Matrix::t(X_b))
			h[idx] <- w[idx] * rowSums(as.matrix(X_b) * t(as.matrix(V_b)))
			if (verbose) utils::setTxtProgressBar(pb, i)
		}
	} else {
		# X^TWX = U^T U (U upper triangular); solve L V = X^T where L = U^T.
		# Then h_i = w_i * ||V[:,i]||^2.
		U <- as(fac, "dtrMatrix")
		V <- forwardsolve(t(as.matrix(U)), t(as.matrix(X))) # p×n
		h <- w * colSums(V^2)
	}

	h[w == 0] <- 0
	h
}

#' @inherit stats::influence
#' @param batch_size integer; passed to `hatvalues.glm4()`.
#' @param verbose logical; passed to `hatvalues.glm4()`.
#' @export
#' @keywords internal
influence.glm4 <- function(model, batch_size = 1000L, verbose = FALSE, ...) {
	h <- hatvalues.glm4(model, batch_size, verbose)
	disp <- dispersion.glm4(model)
	Vmu <- model$family$variance(model$fitted.values)
	pear <- (model$y - model$fitted.values) * sqrt(model$prior.weights) / sqrt(Vmu)
	dev.res <- as.numeric(MatrixModels:::residuals(model$glm4_fit, type = "deviance"))
	list(hat = h, pear.res = pear, dev.res = dev.res, dispersion = disp)
}

# rstandard and cooks.distance broadly follow stats package

#' @inherit stats::rstandard
#' @export
#' @keywords internal
rstandard.glm4 <- function(model,
		infl = influence.glm4(model),
		type = c("pearson", "deviance"), ...) {
	type <- match.arg(type)
	res <- switch(type, pearson = infl$pear.res, infl$dev.res)
	r <- res / sqrt(infl$dispersion * (1 - infl$hat))
	r[is.infinite(r)] <- NaN
	r
}

#' @inherit stats::cooks.distance
#' @export
#' @keywords internal
cooks.distance.glm4 <- function(model, infl = influence.glm4(model), ...) {
	h <- infl$hat
	D <- (infl$pear.res^2 / (infl$dispersion * model$rank)) * (h / (1 - h)^2)
	D[is.infinite(D)] <- NaN
	D
}
