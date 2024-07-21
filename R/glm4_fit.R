#' Fitting Generalized Linear Models Using Sparse Matrices
#'
#' @description
#' `glm4()`, is used to fit generalized linear models,  specified by giving a symbolic description of the linear predictor and a description of the error distribution.
#'
#' It is very similar to the standard `stats::glm()` function, but supports sparse matrices via the Matrix package, which can dramatically improve memory and computational efficiency on large and/or high-dimensional data.
#'
#' @usage glm4(formula, family, data, weights, subset, na.action,
#' start = NULL, etastart, mustart, offset,
#' sparse = FALSE, drop.unused.levels = FALSE, doFit = TRUE,
#' control = list(\dots),
#' model = TRUE, x = FALSE, y = TRUE, contrasts = NULL, \dots)
#'
#' @inheritParams MatrixModels::glm4
#'
#' @details
#' This function is a wrapper for `MatrixModels::glm4()` which returns an more user-friendly object designed to resemble `stats::glm()` as closely as possible. Behind the scenes, it extracts the relevant model details from the S4 class `MatrixModels::glm4()` object and calculates new ones where necessary (e.g. AIC, deviance, residual degrees of freedom, etc) as per `stats::glm()`.
#'
#' @returns A list object of class `glm4`. See `stats::glm()` for more details on returned components.
#' @export

glm4 <- function(formula, data, ...){

	fit <- MatrixModels::glm4(formula, data, ...)

	# Get coef and resid via MM
	coefficients <- MatrixModels::coef(fit)
	residuals <- MatrixModels::residuals(fit, type = "working")

	# Get family functions:
	family <- fit@resp@family
	variance <- family$variance
	linkinv  <- family$linkinv
	dev.resids <- family$dev.resids
	aic <- family$aic
	mu.eta <- family$mu.eta

	# Model outputs
	## respModule
	mu <- fit@resp@mu
	eta <- fit@resp@eta
	weights <- fit@resp@weights
	wt <- c(fit@resp@sqrtXwt)^2
	y <- fit@resp@y
	offset <- fit@resp@offset

	## predModule
	X <- fit@pred@X

	# Model characteristics
	rank <- rank.glm4(fit)
	call <- match.call()
	terms <- stats::terms.formula(formula, data = data)
	intercept <- attr(terms, "intercept") > 0L
	contrasts <- attr(X, "contrasts")

	# Res and null DF calculations (as per glm.fit)
	Xdims <- dim(X)
	nobs <- Xdims[1]
	n.ok <- nobs - sum(weights==0)
	nulldf <- n.ok - as.integer(intercept)
	resdf <- n.ok - rank

	## calculate null deviance -- corrected in glm() if offset and intercept
	wtdmu <- if (intercept)
			sum(weights * y)/sum(weights)
		else
			linkinv(offset)
	nulldev <- sum(dev.resids(y, wtdmu, weights))
	dev <- sum(dev.resids(y, mu, weights))

	# For binomial models, get n, as called in binomial()$initialize
	if (family$family == "binomial"){
		if (NCOL(y) == 1) {
			n <- rep.int(1, nobs)
		}
		else if (NCOL(y) == 2) {
			n <- (y1 <- y[, 1L]) + y[, 2L]
		}
		else stop(
			gettextf(
				"for the '%s' family, y must be a vector of 0 and 1's\nor a 2 column matrix where col 1 is no. successes and col 2 is no. failures",
				"binomial"),
			domain = NA)
	}

	# Calculate AIC
	aic.model <- as.numeric(aic(y, n, mu, weights, dev) + 2*rank)

	# Populate list object and return as class "glm4"
	fit <- list(glm4_fit = fit, coefficients = coefficients, residuals = residuals, fitted.values = mu,
							family = family, linear.predictors = eta, null.deviance = nulldev, deviance = dev,
							iter = fit@fitProps$iter, offset = offset, method = "MatrixModels::glm4",
							contrasts = contrasts, df.residual = resdf, df.null = nulldf, model = X,
							data = data, terms = terms, y = y, call = call, formula = formula,
							aic = aic.model, rank = rank, prior.weights = weights, weights = wt
	)
	class(fit) <- c("glm4")
	return(fit)
}
