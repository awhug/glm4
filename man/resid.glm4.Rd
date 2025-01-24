% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glm4_methods.R
\name{resid.glm4}
\alias{resid.glm4}
\title{Accessing Generalized Linear Model Fits}
\usage{
\method{resid}{glm4}(object, ...)
}
\arguments{
\item{object}{an object of class \code{glm}, typically the result of
    a call to \code{\link[stats]{glm}}.}

\item{...}{further arguments passed to or from other methods.}
}
\description{
These functions are all \code{\link{methods}} for class \code{glm} or
  \code{summary.glm} objects.
}
\details{
The references define the types of residuals: Davison & Snell is a
  good reference for the usages of each.

  The partial residuals are a matrix of working residuals, with each
  column formed by omitting a term from the model.

  How \code{residuals} treats cases with missing values in the original
  fit is determined by the \code{na.action} argument of that fit.
  If \code{na.action = na.omit} omitted cases will not appear in the
  residuals, whereas if \code{na.action = na.exclude} they will appear,
  with residual value \code{NA}.  See also \code{\link[stats]{naresid}}.

  For fits done with \code{y = FALSE} the response values are computed
  from other components.
}
\references{
Davison, A. C. and Snell, E. J. (1991)
  \emph{Residuals and diagnostics.}  In: Statistical Theory
  and Modelling. In Honour of Sir David Cox, FRS, eds.
  Hinkley, D. V., Reid, N. and Snell, E. J., Chapman & Hall.

  Hastie, T. J. and Pregibon, D. (1992)
  \emph{Generalized linear models.}
  Chapter 6 of \emph{Statistical Models in S}
  eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.

  McCullagh P. and Nelder, J. A. (1989)
  \emph{Generalized Linear Models.}
  London: Chapman and Hall.
}
\seealso{
\code{\link[stats]{glm}} for computing \code{glm.obj}, \code{\link[stats]{anova.glm}};
  the corresponding \emph{generic} functions, \code{\link[stats]{summary.glm}},
  \code{\link[stats]{coef}}, \code{\link[stats]{deviance}},
  \code{\link[stats]{df.residual}},
  \code{\link[stats]{effects}}, \code{\link[stats]{fitted}},
  \code{\link[stats]{residuals}}.

  \link[stats]{influence.measures} for deletion diagnostics, including
  standardized (\code{\link[stats]{rstandard}})
  and studentized (\code{\link[stats]{rstudent}}) residuals.
}
\keyword{internal}
