\name{BloodPaTH_model}
\alias{BloodPaTH_model}
\title{
Simulate bloodborne pathogen transmission within a healthcare setting
}
\description{
Agent-based model (coded in C++) that simulates (i) the movements of patients
between wards, (ii) the dynamics of medical devices (contamination, use and
reuse, and sterile-stock availability), and (iii) the dynamics of pathogen
transmission in patients undergoing invasive procedures (following an SEI model).
}
\usage{
BloodPaTH_model(time_step, t, pathogen, min_e_phase, max_e_phase, nb_patients,
  nb_wards, nb_adm, adm_prob, init_prob, prev_type, prev_init, WT_matrix,
  nb_procedures, procedure_names, nb_devices, device_names, nb_devices_new,
  nb_devices_used, nb_devices_contaminated, refill_quantities, refill_freq,
  table_procedures_devices, sterilization_prob, PPM_matrix, dist_risk,
  intervention = "none", prob_screening = 0, at_risk_wards = 0L,
  output = "simple", id_sim = 1L)
}
\arguments{
  \item{time_step}{Time-step, expressed in hours (e.g. 1/60 for a 1-minute time-step). Type: float.}
  \item{t}{Simulation time, expressed in hours (e.g. 24*365 for one year). Type: float.}
  \item{pathogen}{Bloodborne pathogen type: \code{"HCV"}, \code{"HBV"} or \code{"custom"}. If \code{"custom"}, \code{dist_risk} must be informed. \code{"HIV"} is not yet supported. Type: string.}
  \item{min_e_phase}{Minimum value of the eclipse phase, expressed in hours. Type: float.}
  \item{max_e_phase}{Maximum value of the eclipse phase, expressed in hours. Type: float.}
  \item{nb_patients}{Number of patients. Type: integer.}
  \item{nb_wards}{Number of wards within the healthcare setting. Type: integer.}
  \item{nb_adm}{Number of admission routes (typically, the number of departments). Type: integer.}
  \item{adm_prob}{Probability of admission for each admission route; must sum to 1 and be the same length as \code{nb_adm}. Type: float vector.}
  \item{init_prob}{Probability of first admission in each ward, per admission route. Matrix of size \code{nb_adm} x \code{nb_wards}. Type: float matrix.}
  \item{prev_type}{Level at which the initial prevalence is informed: \code{"admission route"} or \code{"ward"}. Type: string.}
  \item{prev_init}{Prevalence upon admission, matching the length of \code{nb_adm} (if \code{prev_type = "admission route"}) or \code{nb_wards} (if \code{prev_type = "ward"}). Type: float vector.}
  \item{WT_matrix}{Merged ward-transition matrices for all admission routes, of size \code{(nb_wards+1)} x \code{(nb_wards+1)} (the extra row/column corresponds to the "discharged" event). See \code{\link{list_to_combined_matrices}}. Type: float matrix.}
  \item{nb_procedures}{Number of procedure types performed within the setting, including the "no procedure" event (e.g. 10 procedures + 1 = 11). Type: integer.}
  \item{procedure_names}{Vector of procedure names. Type: string vector.}
  \item{nb_devices}{Number of different device types used within the setting. Type: integer.}
  \item{device_names}{Vector of device names. Type: string vector.}
  \item{nb_devices_new}{Initial number of new sterile devices in each ward. Matrix of size \code{nb_wards} x \code{nb_devices}. Type: integer matrix.}
  \item{nb_devices_used}{Initial number of non-sterile (previously used) devices in each ward. Matrix of size \code{nb_wards} x \code{nb_devices}. Type: integer matrix.}
  \item{nb_devices_contaminated}{Initial number of non-sterile and contaminated devices in each ward. Matrix of size \code{nb_wards} x \code{nb_devices}. Type: integer matrix.}
  \item{refill_quantities}{Quantity of devices added to the sterile stock at each refill event. Matrix of size \code{nb_wards} x \code{nb_devices}. Type: integer matrix.}
  \item{refill_freq}{Refill frequency for each device type in each ward, expressed in hours. Matrix of size \code{nb_wards} x \code{nb_devices}. Type: integer matrix.}
  \item{table_procedures_devices}{Correspondence table between procedures (column 1) and devices (column 2), one row per procedure/device association. Type: integer matrix.}
  \item{sterilization_prob}{Probability of efficient sterilization for each device type. Type: float vector.}
  \item{PPM_matrix}{Merged matrices of probability of undergoing each procedure (columns) within each ward (rows), for all admission routes. Size \code{nb_wards} x \code{nb_procedures}. See \code{\link{list_to_combined_matrices}}. Type: float matrix.}
  \item{dist_risk}{Data frame giving the distribution parameters of the risk of infection for each procedure. Columns 1-3 hold the distribution parameters and column 4 the distribution name (\code{"lnorm"}, \code{"norm"} or \code{"pert"}); must have \code{nb_procedures} rows, the last row (risk = 0) corresponding to the "no procedure" event. Type: data frame.}
  \item{intervention}{Type of intervention to reduce the burden of infections: \code{"none"}, \code{"upon-admission"} or \code{"ward-level"}. Type: string. Default: \code{"none"}.}
  \item{prob_screening}{Probability of patient screening upon admission, at a setting or ward level. Type: float. Default: \code{0}.}
  \item{at_risk_wards}{Vector of the most at-risk wards (as integer ranks/IDs) to target with a ward-level intervention. Type: integer vector. Default: \code{0}.}
  \item{output}{Either \code{"simple"} or \code{"detailed"} (the latter not yet implemented). Type: string. Default: \code{"simple"}.}
  \item{id_sim}{Simulation ID, used to allow replication and to keep baseline/intervention scenarios comparable across simulations. Type: integer. Default: \code{1}.}
}
\value{
A list (when \code{output = "simple"}) with the following elements:
  \item{incidence}{Number of new cases at each time-step, at a setting level, for the baseline scenario.}
  \item{incidence_intervention}{Same as \code{incidence}, for the intervention scenario (equal to \code{incidence} if \code{intervention = "none"}).}
  \item{susceptible_patients}{Number of susceptible patients at each time-step, at a setting level.}
  \item{susceptible_patients_intervention}{Same as \code{susceptible_patients}, for the intervention scenario.}
  \item{wards_events}{Total number of contamination events in each ward, across the whole period.}
  \item{wards_events_intervention}{Same as \code{wards_events}, for the intervention scenario.}
  \item{wards_susceptibles}{Number of susceptible patients in each ward at each time-step.}
  \item{wards_susceptibles_intervention}{Same as \code{wards_susceptibles}, for the intervention scenario.}
  \item{contaminated_devices}{Total number of contaminated devices of each type at the end of the period.}
  \item{contaminated_devices_intervention}{Same as \code{contaminated_devices}, for the intervention scenario.}
  \item{used_devices}{Total number of non-sterile devices of each type at the end of the period.}
  \item{used_devices_intervention}{Same as \code{used_devices}, for the intervention scenario.}
  \item{infection_events_devices}{Total number of contamination events attributable to each device type.}
  \item{infection_events_devices_intervention}{Same as \code{infection_events_devices}, for the intervention scenario.}
  \item{count_patients_hospital}{Total number of patients that entered the setting across the period.}
  \item{newly_admitted_patients}{Number of newly admitted patients at each time-step.}
  \item{count_admissions_wards}{Number of accesses (not unique) to each ward across the period.}
  \item{count_tests}{Total number of tests performed across the period.}
  \item{newly_contaminated_patients}{Number of newly contaminated patients at each time-step.}
  \item{newly_contaminated_patients_intervention}{Same as \code{newly_contaminated_patients}, for the intervention scenario.}
}
\references{
Henriot P. et al. An agent-based model to simulate the transmission dynamics of
bloodborne pathogens within hospitals. \emph{PLOS Computational Biology}.
\url{https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1012850}
}
\seealso{
\code{\link{list_to_combined_matrices}}, \code{\link{yearly_cum_inc}}, \code{\link{plot_incidence}}
}
\examples{
  \dontrun{
     ## See Application/model_application_example.R on GitHub for a full
     ## worked example on synthetic hospital data.
     output <- BloodPaTH_model(time_step = 1/12, t = 8760, pathogen = "HCV", ...)
  }
}
