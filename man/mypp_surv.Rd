% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aaa.R
\name{mypp_surv}
\alias{mypp_surv}
\title{Find possible parents for each variable based on FDR or correlation
cutoff for phenotype driven two-level possible parents of survival outcome variable.}
\usage{
mypp_surv(mydata, Covariates, alpha1, alpha2)
}
\arguments{
\item{mydata}{Data frame of data with last two colums of survival outcome and time.}

\item{Covariates}{The column names of Covariates in mydata.}

\item{alpha1}{Cutoff value of FDR or correlation for second-level parents for outcome variable.}

\item{y}{Outcome variable.}
}
\value{
a list of possible parents for each variable.
}
\description{
Find possible parents for each variable based on FDR or correlation
cutoff for phenotype driven two-level possible parents of survival outcome variable.
}
