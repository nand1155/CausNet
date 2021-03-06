% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aaa.R
\name{mscores}
\alias{mscores}
\title{Compute the scores for nodes with no parents for each node.}
\usage{
mscores(vset, mydata, surv = 0, scoreFn)
}
\arguments{
\item{vset}{List of Covariates.}

\item{mydata}{Data frame of data with last column as output variable in case of binary/continuous outcome;
in case of survival data, data with last two columns the survival outcome and time.}

\item{surv}{Whether Survival outcome data, 0 or 1.}

\item{scoreFn}{Score Function to be used. Values 'bic' or 'bge'.}
}
\value{
List of scores for nodes with no parents for each node.
}
\description{
Compute the scores for nodes with no parents for each node.
}
