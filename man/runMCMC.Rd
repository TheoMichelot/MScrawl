% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runMCMC.R
\name{runMCMC}
\alias{runMCMC}
\title{Run MCMC iterations}
\usage{
runMCMC(track, nbStates, nbIter, inits, priors, props, tunes, kalmanpars,
  updateState = TRUE)
}
\arguments{
\item{track}{Dataframe of data, with columns "x", "y", and "time"}

\item{nbStates}{Number of states}

\item{nbIter}{Number of iterations}

\item{inits}{List of initial parameters 
(beta, sigma, Q, state)}

\item{priors}{List of parameters of prior distributions, with components:
\itemize{
  \item{"mean":} Vector of means for normal priors on movement parameters, of length 2*nbStates
  \item{"sd":} Vector of standard deviations for normal priors on movement parameters, of length
2*nbStates
  \item{"shape":} Vector of shapes of gamma priors for the transition rates
  \item{"rate":} Vector of rates of gamma priors for the transition rates
  \item{"con":} Vector of concentrations of Dirichlet priors for transition probabilities
}}

\item{props}{List of parameters of proposal distributions, with components:
\itemize{
  \item{"betaSD":} Scalar standard deviation for normal proposal distribution of beta
  \item{"sigmaSD":} Scalar standard deviation for normal proposal distribution of sigma
  \item{"updateLim":} Vector of two values: min and max length of updated state sequence
  \item{"updateProbs":} Probability for each element of updateLim[1]:updateLim[2] (if NULL,
  all values are equiprobable)
}}

\item{tunes}{List of tuning parameters, with components:
\itemize{
  \item{"thinStates":} Thinning factor for the posterior state sequences (needed because
of memory limitations)
}}

\item{kalmanpars}{List of parameters of the Kalman filter, with components:
\itemize{
  \item{"Hmat":} Matrix of observation error variance (four columns, and one row 
for each row of data)
  \item{"a0":} Initial state estimate vector
  \item{"P0":} Initial estimate covariance matrix
}}

\item{updateState}{Logical. If FALSE, the state process is not updated
(for exploratory analysis only)}
}
\description{
Run MCMC iterations
}
