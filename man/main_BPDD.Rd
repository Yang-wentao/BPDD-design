% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_analysis_mainfunction.R
\name{main_BPDD}
\alias{main_BPDD}
\title{Main function for BPDD: A Bayesian phase I/II platform design with data augmentation
accounting for delayed outcomes}
\usage{
main_BPDD(
  n.cohortsize,
  cohortsize,
  phi.T,
  phi.E,
  pt.mat,
  prior.pt.mat,
  prior.pe.mat,
  vpt,
  vd,
  type = "binary",
  C_T = 0.05,
  C_E = 0.05
)
}
\arguments{
\item{n.cohortsize}{The total number of cohorts to simulate, default is 14 or 12.}

\item{cohortsize}{The size of each cohort, default is 3.}

\item{phi.T}{The upper bound of the acceptable toxicity level in a drug trial.}

\item{phi.E}{The lower bound of the acceptable efficacy level in a drug trial.}

\item{pt.mat}{A matrix representing true toxicity values.}

\item{prior.pt.mat}{A matrix representing the priors of toxicity.}

\item{prior.pe.mat}{A matrix representing the priors of efficacy.}

\item{vpt}{A matrix representing the true parameters in setting models.}

\item{vd}{A matrix representing the initial doses from 0 to 1.}

\item{type}{The type outcomes, default is "binary".}

\item{C_T}{The minimum acceptable probability that the drug's toxicity is below the upper bound in the dose-finding process, with a default value of 0.05.}

\item{C_E}{The minimum acceptable probability that the drug's efficacy is above the lower bound in the dose-finding process, with a default value of 0.05.}
}
\value{
A list containing the results of the Bayesian analysis, including list of OBD_dose_pick, parameter_estimate, finish_day, Patient Enrollment Schedule
}
\description{
This function performs the main analysis for Bayesian data augmentation
in the context of delayed outcomes in Phase I/II clinical trials. The
function handles cohort sizes, treatment effects, prior distributions, and
other parameters to model delayed clinical outcomes. It applies Bayesian
methods to simulate and augment data, supporting robust analysis of clinical trial data.
}
