#' @inherit stats::anova.glm
#' @param object an object of class "glm4"
#' @param ... additional objects of class "glm4" for multi-model comparison
#' @export
anova.glm4 <- function(object, ..., dispersion = NULL, test = NULL) {
	dotargs <- list(...)
	named <- if (is.null(names(dotargs)))
		rep_len(FALSE, length(dotargs))
	else (names(dotargs) != "")
	if (any(named))
		warning("the following arguments to 'anova.glm4' are invalid and dropped: ",
						paste(deparse(dotargs[named]), collapse = ", "))
	dotargs <- dotargs[!named]
	is.glm4 <- vapply(dotargs, function(x) inherits(x, "glm4"), NA)
	dotargs <- dotargs[is.glm4]
	if (length(dotargs))
		return(.anova_glm4_multimodel(c(list(object), dotargs), dispersion = dispersion, test = test))

	varlist <- attr(object$terms, "variables")
	response <- as.character(varlist[-1L])[1L]
	tl <- attr(object$terms, "term.labels")
	nvars <- length(tl)
	n <- NROW(object$y)
	resdev <- resdf <- NULL

	if (nvars > 1) {
		y <- object$y
		if (is.null(y)) {
			mu.eta <- object$family$mu.eta
			eta <- object$linear.predictors
			y <- object$fitted.values + object$residuals * mu.eta(eta)
		}

		# Pre-extract dev.resids for use inside the loop
		dev.resids <- object$family$dev.resids

		# Reconstruct a fresh family object from the stored name+link.
		family_ <- do.call(object$family$family, list(link = object$family$link))
		is_sparse <- is(object$model, "sparseMatrix")

		for (i in seq_len(nvars - 1L)) {
			# Refit using a sub-formula on the original data so that
			# MatrixModels builds its own sparse model matrix internally.
			sub_formula <- reformulate(tl[seq_len(i)], response = response)

			# Define fit arguments to pass to MatrixModels::glm4
			fit_args <- list(formula = sub_formula, data = object$data,
							 family = family_, sparse = is_sparse)
			if (length(object$prior.weights) == n)
				fit_args$weights <- object$prior.weights
			if (length(object$offset) == n && any(object$offset != 0))
				fit_args$offset <- object$offset
			fit_s4 <- do.call(MatrixModels::glm4, fit_args)

			# Extract deviance from S4 slots (mirrors glm4_fit)
			resdev <- c(resdev, sum(dev.resids(fit_s4@resp@y, fit_s4@resp@mu, fit_s4@resp@weights)))
			resdf <- c(resdf, df.residual.glm4(fit_s4))
		}
	}

	# Extract residual df and deviance
	resdf <- c(object$df.null, resdf, object$df.residual)
	resdev <- c(object$null.deviance, resdev, object$deviance)

	table <- data.frame(
		Df = c(NA, -diff(resdf)),
		Deviance = c(NA, pmax(0, -diff(resdev))),
		`Resid. Df` = resdf,
		`Resid. Dev` = resdev,
		check.names = FALSE
	)
	if (length(tl) == 0L) table <- table[1, , drop = FALSE]
	rownames(table) <- c("NULL", tl)

	title <- paste0(
		"Analysis of Deviance Table (glm4)",
		"\n\nModel: ", object$family$family,
		", link: ", object$family$link,
		"\n\nResponse: ", response,
		"\n\nTerms added sequentially (first to last)\n\n"
	)

	fam <- object$family

	# Mirror stats::anova.glm default test selection:
	# gaussian (fam$dispersion = NA) -> "F"; binomial/poisson (= 1) -> "Chisq"; else FALSE
	if (is.null(test)) {
		test <- if (!is.null(dispersion))
			"Chisq"
		else if (!is.null(fam$dispersion))
			if (is.na(fam$dispersion)) "F" else "Chisq"
		else FALSE
	}

	if (!isFALSE(test)) {
		if (is.null(dispersion))
			dispersion <- summary(object)$dispersion
		df.dispersion <- if (is.null(fam$dispersion))
			if (isTRUE(dispersion == 1)) Inf else object$df.residual
		else if (is.na(fam$dispersion))
			object$df.residual
		else Inf

		if (isTRUE(test == "F") && df.dispersion == Inf) {
			fname <- fam$family
			if (fname %in% c("binomial", "poisson"))
				warning(gettextf("using F test with a '%s' family is inappropriate",
								 fname), domain = NA)
			else
				warning("using F test with a fixed dispersion is inappropriate")
		}
		table <- stats::stat.anova(
			table = table,
			test = test,
			scale = dispersion,
			df.scale = df.dispersion,
			n = n)
	}
	structure(table, heading = title, class = c("anova", "data.frame"))
}

#' @keywords internal
.anova_glm4_multimodel <- function(object, dispersion = NULL, test = NULL) {
	responses <- vapply(object, function(x) deparse(x$formula[[2L]]), character(1))
	sameresp <- responses == responses[1L]
	if (!all(sameresp)) {
		object <- object[sameresp]
		warning("models with response ", deparse(responses[!sameresp]),
						" removed because response differs from model 1")
	}

	ns <- vapply(object, function(x) NROW(x$model), integer(1))
	if (any(ns != ns[1L]))
		stop("models were not all fitted to the same size of dataset")

	resdf <- vapply(object, `[[`, numeric(1), "df.residual")
	resdev <- vapply(object, `[[`, numeric(1), "deviance")

	table <- data.frame(
		`Resid. Df` = resdf,
		`Resid. Dev` = resdev,
		Df = c(NA, -diff(resdf)),
		Deviance = c(NA, -diff(resdev)),
		check.names = FALSE
	)

	variables <- lapply(object, function(x) paste(deparse(x$formula), collapse = "\n"))
	topnote <- paste("Model", format(seq_along(object)), ":",
										 unlist(variables), sep = " ", collapse = "\n")

	if (!is.null(test)) {
		bigmodel <- object[[length(object)]]
		disp_val <- if (is.null(dispersion)) summary(bigmodel)$dispersion else dispersion
		df.dispersion <- if (disp_val == 1) Inf else resdf[length(resdf)]
		table <- stats::stat.anova(
			table = table, 
			test = test, 
			scale = disp_val, 
			df.scale = df.dispersion, 
			n = ns[1L]
		)
	}

	structure(
		table,
		heading = c("Analysis of Deviance Table (glm4)\n", topnote),
		class = c("anova", "data.frame"))
}
