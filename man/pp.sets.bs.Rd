% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aaa.R
\name{pp.sets.bs}
\alias{pp.sets.bs}
\title{Find the best parent sets and best scores for all sets of possible parents for each node.}
\usage{
pp.sets.bs(pps, ppss, ms, max_parents, surv)
}
\arguments{
\item{pps}{List of Possible parents sets' subsets.}

\item{ppss}{List of Possible parents sets' subsets scores.}

\item{ms}{List of scores of all nodes with no parents.}

\item{max_parents}{Maximum in-degree.}

\item{surv}{Whether Survival outcome data, 0 or 1.}
}
\value{
The list of best parent sets and best scores for all sets of possible parents for each node.
}
\description{
Find the best parent sets and best scores for all sets of possible parents for each node.
}
