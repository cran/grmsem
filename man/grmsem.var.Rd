% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/grmsem.var.R
\name{grmsem.var}
\alias{grmsem.var}
\title{grmsem variance estimation function}
\usage{
grmsem.var(grmsem.out = NULL)
}
\arguments{
\item{grmsem.out}{A grmsem.fit or grmsem.stpar object. Default NULL.}
}
\value{
\code{grmsem.var} returns a list object consisting of the following matrices:
\item{VA}{estimated genetic variance}
\item{VA.se}{standard error of estimated genetic variance}
\item{VE}{estimated residual variance}
\item{VE.se}{standard error of estimated residual variance}
\item{VP}{estimated total phenotypic variance}
\item{RG}{genetic correlation}
\item{RG.se}{standard error genetic correlation}
\item{RE}{residual correlation}
\item{RG.se}{standard error residual correlation}
}
\description{
This function estimates genetic and residual variances, and genetic correlations.
}
\details{
The \code{grmsem.var} function can be used to estimate genetic and residual covariance and correlations for DS, Cholesky, IP and IPC models, based on 
\code{grmsem.fit} or \code{grmsem.stpar} objects. For the latter, the diagonal elements of the VA output matrix 
detail the heritabilities. Except for directly estimated variance components using the DS model, all standard errors 
are derived with the Delta method.
}
\examples{
#(runtime should be less than one minute)
\donttest{
out <- grmsem.fit(ph.small, G.small, LogL = TRUE, estSE = TRUE)
var.out <- grmsem.var(out)
print(var.out)}
}
\keyword{grmsem}
