anova.glm4 <- function (object, ..., dispersion = NULL, test = NULL) {
	dotargs <- list(...)
	named <- if (is.null(names(dotargs)))
		rep_len(FALSE, length(dotargs))
	else (names(dotargs) != "")
	if (any(named))
		warning("the following arguments to 'anova.glm4' are invalid and dropped: ",
						paste(deparse(dotargs[named]), collapse = ", "))
	dotargs <- dotargs[!named]
	is.glm <- vapply(dotargs, function(x) inherits(x, "glm4"),
									 NA)
	dotargs <- dotargs[is.glm]
	if (length(dotargs))
		return(anova.glmlist(c(list(object), dotargs), dispersion = dispersion,
												 test = test))

	varlist <- attr(object$terms, "variables")
	response <- as.character(varlist[-1L])[1L]
	form <- formula(paste(response, "~ . + 0"))
	x <- if (n <- match("x", names(object), 0L))
		object[[n]]
	else model.matrix(object)
	varseq <- attr(x, "assign")
	nvars <- max(0, varseq)
	resdev <- resdf <- NULL

	if (nvars > 1) {
		method <- object$method
		y <- object$y
		if (is.null(y)) {
			mu.eta <- object$family$mu.eta
			eta <- object$linear.predictors
			y <- object$fitted.values + object$residuals * mu.eta(eta)
		} else {
			y_ <- cbind(y)
			colnames(y_) <- response
		}
		for (i in seq_len(max(nvars - 1L, 0))) {
			df_ <- data.frame(y_, x[, varseq <= i, drop = FALSE])
			browser()
			fit <- MatrixModels::glm4(
				formula = form,
				data = df_,
				weights = object$prior.weights,
				start = object$start,
				offset = object$offset,
				family = do.call(object$family$family, list(object$family$link)),
				control = object$control
			)
			resdev <- c(resdev, fit$deviance)
			resdf <- c(resdf, fit$df.residual)
		}
	}
	resdf <- c(object$df.null, resdf, object$df.residual)
	resdev <- c(object$null.deviance, resdev, object$deviance)
	table <- data.frame(c(NA, -diff(resdf)), c(NA, pmax(0, -diff(resdev))),
											resdf, resdev)
	tl <- attr(object$terms, "term.labels")
	if (length(tl) == 0L)
		table <- table[1, , drop = FALSE]
	dimnames(table) <- list(c("NULL", tl), c("Df", "Deviance",
																					 "Resid. Df", "Resid. Dev"))
	title <- paste0("Analysis of Deviance Table", "\n\nModel: ",
									object$family$family, ", link: ", object$family$link,
									"\n\nResponse: ", response, "\n\nTerms added sequentially (first to last)\n\n")
	df.dispersion <- Inf
	if (is.null(dispersion)) {
		dispersion <- summary(object, dispersion = dispersion)$dispersion
		df.dispersion <- if (dispersion == 1)
			Inf
		else object$df.residual
	}
	if (!is.null(test)) {
		if (test == "F" && df.dispersion == Inf) {
			fam <- object$family$family
			if (fam == "binomial" || fam == "poisson")
				warning(gettextf("using F test with a '%s' family is inappropriate",
												 fam), domain = NA)
			else warning("using F test with a fixed dispersion is inappropriate")
		}
		table <- stat.anova(table = table, test = test, scale = dispersion,
												df.scale = df.dispersion, n = NROW(x))
	}
	structure(table, heading = title, class = c("anova", "data.frame"))
}
