\name{BloodPaTH-package}
\alias{BloodPaTH-package}
\alias{BloodPaTH}
\docType{package}
\title{
  Bloodborne Pathogens Transmission in Hospitals
}
\description{
  BloodPaTH is an agent-based model (implemented in C++ via Rcpp/RcppArmadillo) that
  simulates the transmission dynamics of bloodborne pathogens (e.g. HCV, HBV) within
  healthcare settings, using longitudinal retrospective data. It jointly models
  (i) patient movements between wards, (ii) the dynamics of medical devices
  (contamination, use, reuse and sterile-stock availability), and (iii) pathogen
  transmission in patients undergoing invasive procedures, following an SEI model.
}
\details{
  The main entry point is \code{\link{BloodPaTH_model}}, which runs a single
  stochastic simulation and returns the incidence and ward/device-level outputs
  described in its own help page and in the package \href{https://github.com/phenriot/BloodPaTH#readme}{README}.

  Helper functions \code{list_to_combined_matrices}, \code{yearly_cum_inc},
  \code{plot_incidence}, \code{plot_inc_wards} and \code{plot_att_eq} help
  prepare model inputs and summarise/plot simulation outputs across replicates.

  See the \emph{Application} folder on GitHub for a worked example on synthetic data.
}
\author{
Paul Henriot

Maintainer: Paul Henriot <paul.henriot@protonmail.com>
}
\references{
  Henriot P. et al. An agent-based model to simulate the transmission dynamics of
  bloodborne pathogens within hospitals. \emph{PLOS Computational Biology}.
  \url{https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1012850}
}
\keyword{ package }
\seealso{
  \code{\link{BloodPaTH_model}}
}
\examples{
  \dontrun{
     ## See the model_application_example.R script in the Application folder
     ## for a full worked example on synthetic hospital data.
     output <- BloodPaTH_model(time_step = 1/12, t = 8760, pathogen = "HCV", ...)
  }
}
