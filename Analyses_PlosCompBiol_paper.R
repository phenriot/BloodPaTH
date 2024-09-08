########################
### Loading packages ###
########################

library(parallel)
library(graphics)
library(foreach)
library(doParallel)
library(BloodPaTH)

##############################
### Setting up the cluster ###
#############################

n_cores = detectCores()- 6

clust <- makeCluster(n_cores) 

registerDoParallel(clust)

n_sim = 100

##################
### Parameters ###
##################

n_patients = 970
time=105120 # a year in 5-minutes time-steps

Qmax_mat = yearly_eq_use 
N = Qmax_mat/(365*12*24)

## High-resource setting ##
Qn_High = available_eq
R_mat_High = Qn_High/Qmax_mat
R_mat_High[which(is.nan(R_mat_High))]=0
ts_High=round((Qn_High*(R_mat_High-1)+Qmax_mat)/N) ; ts_High[which(is.nan(ts_High))]=0
Q_mat_new_init_High = t(round(R_mat_High*Qn_High))
Q_mat_used_init_High = t(round((1-R_mat_High)*Qn_High)) ; Q_mat_used_init_High[which(Q_mat_used_init_High<0)] =0;
Q_mat_cont_init_High = Q_mat_used_init_High*0

refill_quantities_High=round(t(Qn_High)) # Quantity of device to add to the pool of sterile device for each type at each refill event
refill_freq_High=t(ts_High)#Refill frequency for each device type (must be expressed in hours)

#Sterilization probabilities for each type of device
good_ster_High= rep(0.95,10)

## Low-resource setting ##
Qn_Low = available_eq*0.5
R_mat_Low = (Qn_Low)/Qmax_mat
R_mat_Low[which(is.nan(R_mat_Low))]=0
ts_Low=round((Qn_Low*(R_mat_Low-1)+Qmax_mat)/N) ; ts_Low[which(is.nan(ts_Low))]=0 
Q_mat_new_init_Low = t(round(R_mat_Low*Qn_Low))
Q_mat_used_init_Low = t(round((1-R_mat_Low)*Qn_Low)) ; Q_mat_used_init_Low[which(Q_mat_used_init_Low<0)] =0;
Q_mat_cont_init_Low = Q_mat_used_init_Low*0

refill_quantitie_Low=round(t(Qn_Low)) # Quantity of device to add to the pool of sterile device for each type at each refill event
refill_freq_Low=t(ts_Low)# Refill frequency for each device type (must be expressed in hours)

#Sterilization probabilities for each type of device
good_ster_Low= rep(0.8,10)

################
### HCV Case ###
################

##########################
### Baseline scenarios ###
##########################

min_e_phase_HCV = 576
max_e_phase_HCV = 4032

prev_init = prev_init_ward
prev_type = "ward"
intervention = "none"
at_risk_wards = c(0)

##################################
### High-resource setting case ###
##################################

output_baseline_High = foreach(n = 1:n_sim ,.packages = c("Rcpp","BloodPaTH"))%dopar%{
  
  BloodPaTH_model(t=time, 
                 n_patients=n_patients, 
                 nb_wards=28, 
                 nb_adm=2, 
                 prev_init = prev_init, 
                 prev_type = prev_type,
                 WT_matrix=entry_conversion(list_TM,nb_wards=n_wards,type=1),
                 init_prob = rbind(freq.adm.SH,freq.adm.IM), 
                 adm_prob=c(0.5,0.5),
                 nb_procedures=nb_proc, 
                 nb_equipments=10,
                 PPM_matrix=entry_conversion(list_PPM,nb_procedures=nb_proc,type=3),
                 dist_risk=dist_risk,
                 nb_material_new=Q_mat_new_init_High,
                 nb_material_used=Q_mat_used_init_High,
                 nb_material_contaminated=Q_mat_cont_init_High,
                 min_e_phase=min_e_phase_HCV,
                 max_e_phase=max_e_phase_HCV,
                 time_step=time_step,
                 patient_followup=F,
                 table_proc_equip=as.matrix(corresp_eq_proc),
                 procedure_names=colnames(matrix_proc_eq[,4:18]),
                 equipment_names=matrix_proc_eq$item,
                 sterilization_prob=good_ster_High,
                 eq_quantity = "fixed",
                 refill_quantities=refill_quantities_High,
                 refill_freq=refill_freq_High,
                 equip_bin=matrix(data = 0),
                 output = "simple",
                 prob_screening=0,
                 intervention = intervention,
                 at_risk_wards = at_risk_wards,
                 threshold_reuse= rep(-1,10),
                 id_sim = n)
}


#################################
### Low resource setting case ###
#################################

output_baseline_Low = foreach(n = 1:n_sim ,.packages = c("Rcpp","BloodPaTH"))%dopar%{
  
  set.seed(n)
  
  BloodPaTH_model(t=time, #simulation time
                 n_patients=n_patients, # number of patients
                 nb_wards=28, # number of wards
                 nb_adm=2, # number of admission routes
                 prev_init = prev_init, # initial prevalence for each admission route
                 prev_type = prev_type,
                 WT_matrix=entry_conversion(list_TM,nb_wards=n_wards,type=1), #transition matrix, specific format (merged matrix for all admission routes) but any time-step will do
                 init_prob = rbind(freq.adm.SH,freq.adm.IM), # matrix of probabilities of admission for first ward for both adm routes, nb_adm rows . 
                 adm_prob=c(0.5,0.5),# probabilities of admission for each route
                 nb_procedures=nb_proc, #number of procedures 
                 nb_equipments=10, # number of equipments
                 PPM_matrix=entry_conversion(list_PPM,nb_procedures=nb_proc,type=3), # probabilities of undergoing procedures (columns) in each ward (row), merged for both admission routes; Any time-step is ok but need to be the same as the transition matrix
                 dist_risk=dist_risk, # distribution parameters of the risk for each procedure (dist can be normal, log-normal or pert)
                 nb_material_new=Q_mat_new_init_Low , # quantity of equipment for each type in  each ward
                 nb_material_used=Q_mat_used_init_Low, # quantity of equipment for each type in  each ward
                 nb_material_contaminated=Q_mat_cont_init_Low, # quantity of equipment for each type in  each ward
                 min_e_phase=min_e_phase_HCV, # min value for the eclipse phase
                 max_e_phase=max_e_phase_HCV, # max value for the eclipse phase
                 time_step=time_step, # must be expressed in hours (ex: if time step is 5min then inform 1/12)
                 patient_followup=F, # useless
                 table_proc_equip=as.matrix(corresp_eq_proc), # correspondence table between procedures (col 1) and equipments (col 2)
                 procedure_names=colnames(matrix_proc_eq[,4:18]), # vector of procedure names 
                 equipment_names=matrix_proc_eq$item, # vector of equipments names
                 sterilization_prob=good_ster_Low, # probability of good sterilization for each equipment
                 eq_quantity = "fixed",
                 refill_quantities=refill_quantities_Low, # quantity of equipment to add to the pool of unused equipment for each type at each refill event
                 refill_freq=refill_freq_Low,# refill frequency for each equipment type (must be expressed in hours)
                 equip_bin=round(matrix(nrow = 10,ncol = n_wards,data = runif(280)*100000)),
                 output = "simple",
                 prob_screening=0,
                 intervention = intervention,
                 at_risk_wards = at_risk_wards,
                 threshold_reuse= rep(-1,10),
                 pathogen = "HCV",
                 id_sim = n)
}

################################
### Ward-based interventions ###
################################

prev_init = prev_init_ward
prev_type = "ward"
intervention = "ward-based"

#######################################
### High-resource setting case case ###
#######################################

at_risk_wards = c(28,22,13)

output_WB_High = foreach(n = 1:n_sim ,.packages = c("Rcpp","BloodPaTH"))%dopar%{
  
  set.seed(n)
  
  BloodPaTH_model(t=time, #simulation time
                 n_patients=n_patients, # number of patients
                 nb_wards=28, # number of wards
                 nb_adm=2, # number of admission routes
                 prev_init = prev_init, # initial prevalence for each admission route
                 prev_type = prev_type,
                 WT_matrix=entry_conversion(list_TM,nb_wards=n_wards,type=1), #transition matrix, specific format (merged matrix for all admission routes) but any time-step will do
                 init_prob = rbind(freq.adm.SH,freq.adm.IM), # matrix of probabilities of admission for first ward for both adm routes, nb_adm rows . 
                 adm_prob=c(0.5,0.5),# probabilities of admission for each route
                 nb_procedures=nb_proc, #number of procedures 
                 nb_equipments=10, # number of equipments
                 PPM_matrix=entry_conversion(list_PPM,nb_procedures=nb_proc,type=3), # probabilities of undergoing procedures (columns) in each ward (row), merged for both admission routes; Any time-step is ok but need to be the same as the transition matrix
                 dist_risk=dist_risk, # distribution parameters of the risk for each procedure (dist can be normal, log-normal or pert)
                 nb_material_new=Q_mat_new_init_High , # quantity of equipment for each type in  each ward
                 nb_material_used=Q_mat_used_init_High, # quantity of equipment for each type in  each ward
                 nb_material_contaminated=Q_mat_cont_init_High, # quantity of equipment for each type in  each ward
                 min_e_phase=min_e_phase_HCV, # min value for the eclipse phase
                 max_e_phase=max_e_phase_HCV, # max value for the eclipse phase
                 time_step=time_step, # must be expressed in hours (ex: if time step is 5min then inform 1/12)
                 patient_followup=F, # useless
                 table_proc_equip=as.matrix(corresp_eq_proc), # correspondence table between procedures (col 1) and equipments (col 2)
                 procedure_names=colnames(matrix_proc_eq[,4:18]), # vector of procedure names 
                 equipment_names=matrix_proc_eq$item, # vector of equipments names
                 sterilization_prob=good_ster_High, # probability of good sterilization for each equipment
                 eq_quantity = "fixed",
                 refill_quantities=refill_quantities_High, # quantity of equipment to add to the pool of unused equipment for each type at each refill event
                 refill_freq=refill_freq_High,# refill frequency for each equipment type (must be expressed in hours)
                 equip_bin=matrix(data = 0),
                 output = "simple",
                 prob_screening=1,
                 intervention = intervention,
                 at_risk_wards = at_risk_wards,
                 threshold_reuse= rep(-1,10),
                 id_sim = n)
}

#################################
### Low-resource setting case ###
#################################

at_risk_wards = c(28,27,22)

output_WB_Low = foreach(n = 1:n_sim ,.packages = c("Rcpp","BloodPaTH"))%dopar%{
  
  set.seed(n)
  
  BloodPaTH_model(t=time, #simulation time
                 n_patients=n_patients, # number of patients
                 nb_wards=28, # number of wards
                 nb_adm=2, # number of admission routes
                 prev_init = prev_init, # initial prevalence for each admission route
                 prev_type = prev_type,
                 WT_matrix=entry_conversion(list_TM,nb_wards=n_wards,type=1), #transition matrix, specific format (merged matrix for all admission routes) but any time-step will do
                 init_prob = rbind(freq.adm.SH,freq.adm.IM), # matrix of probabilities of admission for first ward for both adm routes, nb_adm rows . 
                 adm_prob=c(0.5,0.5),# probabilities of admission for each route
                 nb_procedures=nb_proc, #number of procedures 
                 nb_equipments=10, # number of equipments
                 PPM_matrix=entry_conversion(list_PPM,nb_procedures=nb_proc,type=3), # probabilities of undergoing procedures (columns) in each ward (row), merged for both admission routes; Any time-step is ok but need to be the same as the transition matrix
                 dist_risk=dist_risk, # distribution parameters of the risk for each procedure (dist can be normal, log-normal or pert)
                 nb_material_new=Q_mat_new_init_Low , # quantity of equipment for each type in  each ward
                 nb_material_used=Q_mat_used_init_Low, # quantity of equipment for each type in  each ward
                 nb_material_contaminated=Q_mat_cont_init_Low, # quantity of equipment for each type in  each ward
                 min_e_phase=min_e_phase_HCV, # min value for the eclipse phase
                 max_e_phase=max_e_phase_HCV, # max value for the eclipse phase
                 time_step=time_step, # must be expressed in hours (ex: if time step is 5min then inform 1/12)
                 patient_followup=F, # useless
                 table_proc_equip=as.matrix(corresp_eq_proc), # correspondence table between procedures (col 1) and equipments (col 2)
                 procedure_names=colnames(matrix_proc_eq[,4:18]), # vector of procedure names 
                 equipment_names=matrix_proc_eq$item, # vector of equipments names
                 sterilization_prob=good_ster_Low, # probability of good sterilization for each equipment
                 eq_quantity = "fixed",
                 refill_quantities=refill_quantities_Low, # quantity of equipment to add to the pool of unused equipment for each type at each refill event
                 refill_freq=refill_freq_Low,# refill frequency for each equipment type (must be expressed in hours)
                 equip_bin=round(matrix(nrow = 10,ncol = n_wards,data = runif(280)*100000)),
                 output = "simple",
                 prob_screening=1,
                 intervention = intervention,
                 at_risk_wards = at_risk_wards,
                 threshold_reuse= rep(-1,10),
                 pathogen = "HCV",
                 id_sim = n)
}


###################################
### Patient-based interventions ###
###################################

prev_init = prev_init_ward
prev_type = "ward"
intervention = "patient-based"

#######################################
### High-resource setting case case ###
#######################################

at_risk_wards = c(28,27,13)

output_PB_High = foreach(n = 1:n_sim ,.packages = c("Rcpp","BloodPaTH"))%dopar%{
  
  BloodPaTH_model(t=time, #simulation time
                 n_patients=n_patients, # number of patients
                 nb_wards=28, # number of wards
                 nb_adm=2, # number of admission routes
                 prev_init = prev_init, # initial prevalence for each admission route
                 prev_type = prev_type,
                 WT_matrix=entry_conversion(list_TM,nb_wards=n_wards,type=1), #transition matrix, specific format (merged matrix for all admission routes) but any time-step will do
                 init_prob = rbind(freq.adm.SH,freq.adm.IM), # matrix of probabilities of admission for first ward for both adm routes, nb_adm rows . 
                 adm_prob=c(0.5,0.5),# probabilities of admission for each route
                 nb_procedures=nb_proc, #number of procedures 
                 nb_equipments=10, # number of equipments
                 PPM_matrix=entry_conversion(list_PPM,nb_procedures=nb_proc,type=3), # probabilities of undergoing procedures (columns) in each ward (row), merged for both admission routes; Any time-step is ok but need to be the same as the transition matrix
                 dist_risk=dist_risk, # distribution parameters of the risk for each procedure (dist can be normal, log-normal or pert)
                 nb_material_new=Q_mat_new_init_High , # quantity of equipment for each type in  each ward
                 nb_material_used=Q_mat_used_init_High, # quantity of equipment for each type in  each ward
                 nb_material_contaminated=Q_mat_cont_init_High, # quantity of equipment for each type in  each ward
                 min_e_phase=min_e_phase_HCV, # min value for the eclipse phase
                 max_e_phase=max_e_phase_HCV, # max value for the eclipse phase
                 time_step=time_step, # must be expressed in hours (ex: if time step is 5min then inform 1/12)
                 patient_followup=F, # useless
                 table_proc_equip=as.matrix(corresp_eq_proc), # correspondence table between procedures (col 1) and equipments (col 2)
                 procedure_names=colnames(matrix_proc_eq[,4:18]), # vector of procedure names 
                 equipment_names=matrix_proc_eq$item, # vector of equipments names
                 sterilization_prob=good_ster_High, # probability of good sterilization for each equipment
                 eq_quantity = "fixed",
                 refill_quantities=refill_quantities_High, # quantity of equipment to add to the pool of unused equipment for each type at each refill event
                 refill_freq=refill_freq_High,# refill frequency for each equipment type (must be expressed in hours)
                 equip_bin=matrix(data = 0),
                 output = "simple",
                 prob_screening=0.405,
                 intervention = intervention,
                 at_risk_wards = at_risk_wards,
                 threshold_reuse= rep(-1,10),
                 pathogen="HCV",
                 id_sim = n)
}


#################################
### Low-resource setting case ###
#################################

at_risk_wards = c(28,27,22)

output_PB_Low = foreach(n = 1:n_sim ,.packages = c("Rcpp","BloodPaTH"))%dopar%{
  
  BloodPaTH_model(t=time, #simulation time
                 n_patients=n_patients, # number of patients
                 nb_wards=28, # number of wards
                 nb_adm=2, # number of admission routes
                 prev_init = prev_init, # initial prevalence for each admission route
                 prev_type = prev_type,
                 WT_matrix=entry_conversion(list_TM,nb_wards=n_wards,type=1), #transition matrix, specific format (merged matrix for all admission routes) but any time-step will do
                 init_prob = rbind(freq.adm.SH,freq.adm.IM), # matrix of probabilities of admission for first ward for both adm routes, nb_adm rows . 
                 adm_prob=c(0.5,0.5),# probabilities of admission for each route
                 nb_procedures=nb_proc, #number of procedures 
                 nb_equipments=10, # number of equipments
                 PPM_matrix=entry_conversion(list_PPM,nb_procedures=nb_proc,type=3), # probabilities of undergoing procedures (columns) in each ward (row), merged for both admission routes; Any time-step is ok but need to be the same as the transition matrix
                 dist_risk=dist_risk, # distribution parameters of the risk for each procedure (dist can be normal, log-normal or pert)
                 nb_material_new=Q_mat_new_init_Low , # quantity of equipment for each type in  each ward
                 nb_material_used=Q_mat_used_init_Low, # quantity of equipment for each type in  each ward
                 nb_material_contaminated=Q_mat_cont_init_Low, # quantity of equipment for each type in  each ward
                 min_e_phase=min_e_phase_HCV, # min value for the eclipse phase
                 max_e_phase=max_e_phase_HCV, # max value for the eclipse phase
                 time_step=time_step, # must be expressed in hours (ex: if time step is 5min then inform 1/12)
                 patient_followup=F, # useless
                 table_proc_equip=as.matrix(corresp_eq_proc), # correspondence table between procedures (col 1) and equipments (col 2)
                 procedure_names=colnames(matrix_proc_eq[,4:18]), # vector of procedure names 
                 equipment_names=matrix_proc_eq$item, # vector of equipments names
                 sterilization_prob=good_ster_Low, # probability of good sterilization for each equipment
                 eq_quantity = "fixed",
                 refill_quantities=refill_quantities_Low, # quantity of equipment to add to the pool of unused equipment for each type at each refill event
                 refill_freq=refill_freq_Low,# refill frequency for each equipment type (must be expressed in hours)
                 equip_bin=round(matrix(nrow = 10,ncol = n_wards,data = runif(280)*100000)),
                 output = "simple",
                 prob_screening=0.72,
                 intervention = intervention,
                 at_risk_wards = at_risk_wards,
                 threshold_reuse= rep(-1,10),
                 pathogen = "HCV",
                 id_sim = n)
}


################
### HBV Case ###
################

##########################
### Baseline scenarios ###
##########################

min_e_phase_HBV = 4320
max_e_phase_HBV = 4320

prev_init = prev_init_ward_HBV
prev_type = "ward"
intervention = "none"
at_risk_wards = c(0)


##################################
### High-resource setting case ###
##################################

output_baseline_High = foreach(n = 1:n_sim ,.packages = c("Rcpp","BloodPaTH"))%dopar%{
  
  BloodPaTH_model(t=time, 
                 n_patients=n_patients, 
                 nb_wards=28, 
                 nb_adm=2, 
                 prev_init = prev_init, 
                 prev_type = prev_type,
                 WT_matrix=entry_conversion(list_TM,nb_wards=n_wards,type=1),
                 init_prob = rbind(freq.adm.SH,freq.adm.IM), 
                 adm_prob=c(0.5,0.5),
                 nb_procedures=nb_proc, 
                 nb_equipments=10,
                 PPM_matrix=entry_conversion(list_PPM,nb_procedures=nb_proc,type=3),
                 dist_risk=dist_risk,
                 nb_material_new=Q_mat_new_init_High,
                 nb_material_used=Q_mat_used_init_High,
                 nb_material_contaminated=Q_mat_cont_init_High,
                 min_e_phase=min_e_phase_HBV,
                 max_e_phase=max_e_phase_HBV,
                 time_step=time_step,
                 patient_followup=F,
                 table_proc_equip=as.matrix(corresp_eq_proc),
                 procedure_names=colnames(matrix_proc_eq[,4:18]),
                 equipment_names=matrix_proc_eq$item,
                 sterilization_prob=good_ster_High,
                 eq_quantity = "fixed",
                 refill_quantities=refill_quantities_High,
                 refill_freq=refill_freq_High,
                 equip_bin=matrix(data = 0),
                 output = "simple",
                 prob_screening=0,
                 intervention = intervention,
                 at_risk_wards = at_risk_wards,
                 threshold_reuse= rep(-1,10),
                 id_sim = n)
}


#################################
### Low resource setting case ###
#################################

output_baseline_Low = foreach(n = 1:n_sim ,.packages = c("Rcpp","BloodPaTH"))%dopar%{
  
  set.seed(n)
  
  BloodPaTH_model(t=time, #simulation time
                 n_patients=n_patients, # number of patients
                 nb_wards=28, # number of wards
                 nb_adm=2, # number of admission routes
                 prev_init = prev_init, # initial prevalence for each admission route
                 prev_type = prev_type,
                 WT_matrix=entry_conversion(list_TM,nb_wards=n_wards,type=1), #transition matrix, specific format (merged matrix for all admission routes) but any time-step will do
                 init_prob = rbind(freq.adm.SH,freq.adm.IM), # matrix of probabilities of admission for first ward for both adm routes, nb_adm rows . 
                 adm_prob=c(0.5,0.5),# probabilities of admission for each route
                 nb_procedures=nb_proc, #number of procedures 
                 nb_equipments=10, # number of equipments
                 PPM_matrix=entry_conversion(list_PPM,nb_procedures=nb_proc,type=3), # probabilities of undergoing procedures (columns) in each ward (row), merged for both admission routes; Any time-step is ok but need to be the same as the transition matrix
                 dist_risk=dist_risk, # distribution parameters of the risk for each procedure (dist can be normal, log-normal or pert)
                 nb_material_new=Q_mat_new_init_Low , # quantity of equipment for each type in  each ward
                 nb_material_used=Q_mat_used_init_Low, # quantity of equipment for each type in  each ward
                 nb_material_contaminated=Q_mat_cont_init_Low, # quantity of equipment for each type in  each ward
                 min_e_phase=min_e_phase_HBV, # min value for the eclipse phase
                 max_e_phase=max_e_phase_HBV, # max value for the eclipse phase
                 time_step=time_step, # must be expressed in hours (ex: if time step is 5min then inform 1/12)
                 patient_followup=F, # useless
                 table_proc_equip=as.matrix(corresp_eq_proc), # correspondence table between procedures (col 1) and equipments (col 2)
                 procedure_names=colnames(matrix_proc_eq[,4:18]), # vector of procedure names 
                 equipment_names=matrix_proc_eq$item, # vector of equipments names
                 sterilization_prob=good_ster_Low, # probability of good sterilization for each equipment
                 eq_quantity = "fixed",
                 refill_quantities=refill_quantities_Low, # quantity of equipment to add to the pool of unused equipment for each type at each refill event
                 refill_freq=refill_freq_Low,# refill frequency for each equipment type (must be expressed in hours)
                 equip_bin=round(matrix(nrow = 10,ncol = n_wards,data = runif(280)*100000)),
                 output = "simple",
                 prob_screening=0,
                 intervention = intervention,
                 at_risk_wards = at_risk_wards,
                 threshold_reuse= rep(-1,10),
                 pathogen = "HBV",
                 id_sim = n)
}

################################
### Ward-based interventions ###
################################

prev_init = prev_init_ward
prev_type = "ward"
intervention = "ward-based"

#######################################
### High-resource setting case case ###
#######################################

at_risk_wards = c(28,22,13)

output_WB_High = foreach(n = 1:n_sim ,.packages = c("Rcpp","BloodPaTH"))%dopar%{
  
  set.seed(n)
  
  BloodPaTH_model(t=time, #simulation time
                 n_patients=n_patients, # number of patients
                 nb_wards=28, # number of wards
                 nb_adm=2, # number of admission routes
                 prev_init = prev_init, # initial prevalence for each admission route
                 prev_type = prev_type,
                 WT_matrix=entry_conversion(list_TM,nb_wards=n_wards,type=1), #transition matrix, specific format (merged matrix for all admission routes) but any time-step will do
                 init_prob = rbind(freq.adm.SH,freq.adm.IM), # matrix of probabilities of admission for first ward for both adm routes, nb_adm rows . 
                 adm_prob=c(0.5,0.5),# probabilities of admission for each route
                 nb_procedures=nb_proc, #number of procedures 
                 nb_equipments=10, # number of equipments
                 PPM_matrix=entry_conversion(list_PPM,nb_procedures=nb_proc,type=3), # probabilities of undergoing procedures (columns) in each ward (row), merged for both admission routes; Any time-step is ok but need to be the same as the transition matrix
                 dist_risk=dist_risk, # distribution parameters of the risk for each procedure (dist can be normal, log-normal or pert)
                 nb_material_new=Q_mat_new_init_High , # quantity of equipment for each type in  each ward
                 nb_material_used=Q_mat_used_init_High, # quantity of equipment for each type in  each ward
                 nb_material_contaminated=Q_mat_cont_init_High, # quantity of equipment for each type in  each ward
                 min_e_phase=min_e_phase_HBV, # min value for the eclipse phase
                 max_e_phase=max_e_phase_HBV, # max value for the eclipse phase
                 time_step=time_step, # must be expressed in hours (ex: if time step is 5min then inform 1/12)
                 patient_followup=F, # useless
                 table_proc_equip=as.matrix(corresp_eq_proc), # correspondence table between procedures (col 1) and equipments (col 2)
                 procedure_names=colnames(matrix_proc_eq[,4:18]), # vector of procedure names 
                 equipment_names=matrix_proc_eq$item, # vector of equipments names
                 sterilization_prob=good_ster_High, # probability of good sterilization for each equipment
                 eq_quantity = "fixed",
                 refill_quantities=refill_quantities_High, # quantity of equipment to add to the pool of unused equipment for each type at each refill event
                 refill_freq=refill_freq_High,# refill frequency for each equipment type (must be expressed in hours)
                 equip_bin=matrix(data = 0),
                 output = "simple",
                 prob_screening=1,
                 intervention = intervention,
                 at_risk_wards = at_risk_wards,
                 pathogen="HBV",
                 threshold_reuse= rep(-1,10),
                 id_sim = n)
}

#################################
### Low-resource setting case ###
#################################

at_risk_wards = c(28,27,22)

output_WB_Low = foreach(n = 1:n_sim ,.packages = c("Rcpp","BloodPaTH"))%dopar%{
  
  set.seed(n)
  
  BloodPaTH_model(t=time, #simulation time
                 n_patients=n_patients, # number of patients
                 nb_wards=28, # number of wards
                 nb_adm=2, # number of admission routes
                 prev_init = prev_init, # initial prevalence for each admission route
                 prev_type = prev_type,
                 WT_matrix=entry_conversion(list_TM,nb_wards=n_wards,type=1), #transition matrix, specific format (merged matrix for all admission routes) but any time-step will do
                 init_prob = rbind(freq.adm.SH,freq.adm.IM), # matrix of probabilities of admission for first ward for both adm routes, nb_adm rows . 
                 adm_prob=c(0.5,0.5),# probabilities of admission for each route
                 nb_procedures=nb_proc, #number of procedures 
                 nb_equipments=10, # number of equipments
                 PPM_matrix=entry_conversion(list_PPM,nb_procedures=nb_proc,type=3), # probabilities of undergoing procedures (columns) in each ward (row), merged for both admission routes; Any time-step is ok but need to be the same as the transition matrix
                 dist_risk=dist_risk, # distribution parameters of the risk for each procedure (dist can be normal, log-normal or pert)
                 nb_material_new=Q_mat_new_init_Low , # quantity of equipment for each type in  each ward
                 nb_material_used=Q_mat_used_init_Low, # quantity of equipment for each type in  each ward
                 nb_material_contaminated=Q_mat_cont_init_Low, # quantity of equipment for each type in  each ward
                 min_e_phase=min_e_phase_HBV, # min value for the eclipse phase
                 max_e_phase=max_e_phase_HBV, # max value for the eclipse phase
                 time_step=time_step, # must be expressed in hours (ex: if time step is 5min then inform 1/12)
                 patient_followup=F, # useless
                 table_proc_equip=as.matrix(corresp_eq_proc), # correspondence table between procedures (col 1) and equipments (col 2)
                 procedure_names=colnames(matrix_proc_eq[,4:18]), # vector of procedure names 
                 equipment_names=matrix_proc_eq$item, # vector of equipments names
                 sterilization_prob=good_ster_Low, # probability of good sterilization for each equipment
                 eq_quantity = "fixed",
                 refill_quantities=refill_quantities_Low, # quantity of equipment to add to the pool of unused equipment for each type at each refill event
                 refill_freq=refill_freq_Low,# refill frequency for each equipment type (must be expressed in hours)
                 equip_bin=round(matrix(nrow = 10,ncol = n_wards,data = runif(280)*100000)),
                 output = "simple",
                 prob_screening=1,
                 intervention = intervention,
                 at_risk_wards = at_risk_wards,
                 threshold_reuse= rep(-1,10),
                 pathogen = "HBV",
                 id_sim = n)
}


###################################
### Patient-based interventions ###
###################################

prev_init = prev_init_ward
prev_type = "ward"
intervention = "patient-based"

#######################################
### High-resource setting case case ###
#######################################

at_risk_wards = c(28,27,13)

output_PB_High = foreach(n = 1:n_sim ,.packages = c("Rcpp","BloodPaTH"))%dopar%{
  
  BloodPaTH_model(t=time, #simulation time
                 n_patients=n_patients, # number of patients
                 nb_wards=28, # number of wards
                 nb_adm=2, # number of admission routes
                 prev_init = prev_init, # initial prevalence for each admission route
                 prev_type = prev_type,
                 WT_matrix=entry_conversion(list_TM,nb_wards=n_wards,type=1), #transition matrix, specific format (merged matrix for all admission routes) but any time-step will do
                 init_prob = rbind(freq.adm.SH,freq.adm.IM), # matrix of probabilities of admission for first ward for both adm routes, nb_adm rows . 
                 adm_prob=c(0.5,0.5),# probabilities of admission for each route
                 nb_procedures=nb_proc, #number of procedures 
                 nb_equipments=10, # number of equipments
                 PPM_matrix=entry_conversion(list_PPM,nb_procedures=nb_proc,type=3), # probabilities of undergoing procedures (columns) in each ward (row), merged for both admission routes; Any time-step is ok but need to be the same as the transition matrix
                 dist_risk=dist_risk, # distribution parameters of the risk for each procedure (dist can be normal, log-normal or pert)
                 nb_material_new=Q_mat_new_init_High , # quantity of equipment for each type in  each ward
                 nb_material_used=Q_mat_used_init_High, # quantity of equipment for each type in  each ward
                 nb_material_contaminated=Q_mat_cont_init_High, # quantity of equipment for each type in  each ward
                 min_e_phase=min_e_phase_HBV, # min value for the eclipse phase
                 max_e_phase=max_e_phase_HBV, # max value for the eclipse phase
                 time_step=time_step, # must be expressed in hours (ex: if time step is 5min then inform 1/12)
                 patient_followup=F, # useless
                 table_proc_equip=as.matrix(corresp_eq_proc), # correspondence table between procedures (col 1) and equipments (col 2)
                 procedure_names=colnames(matrix_proc_eq[,4:18]), # vector of procedure names 
                 equipment_names=matrix_proc_eq$item, # vector of equipments names
                 sterilization_prob=good_ster_High, # probability of good sterilization for each equipment
                 eq_quantity = "fixed",
                 refill_quantities=refill_quantities_High, # quantity of equipment to add to the pool of unused equipment for each type at each refill event
                 refill_freq=refill_freq_High,# refill frequency for each equipment type (must be expressed in hours)
                 equip_bin=matrix(data = 0),
                 output = "simple",
                 prob_screening=0.405,
                 intervention = intervention,
                 at_risk_wards = at_risk_wards,
                 threshold_reuse= rep(-1,10),
                 pathogen="HBV",
                 id_sim = n)
}


#################################
### Low-resource setting case ###
#################################

at_risk_wards = c(28,27,22)

output_PB_Low = foreach(n = 1:n_sim ,.packages = c("Rcpp","BloodPaTH"))%dopar%{
  
  BloodPaTH_model(t=time, #simulation time
                 n_patients=n_patients, # number of patients
                 nb_wards=28, # number of wards
                 nb_adm=2, # number of admission routes
                 prev_init = prev_init, # initial prevalence for each admission route
                 prev_type = prev_type,
                 WT_matrix=entry_conversion(list_TM,nb_wards=n_wards,type=1), #transition matrix, specific format (merged matrix for all admission routes) but any time-step will do
                 init_prob = rbind(freq.adm.SH,freq.adm.IM), # matrix of probabilities of admission for first ward for both adm routes, nb_adm rows . 
                 adm_prob=c(0.5,0.5),# probabilities of admission for each route
                 nb_procedures=nb_proc, #number of procedures 
                 nb_equipments=10, # number of equipments
                 PPM_matrix=entry_conversion(list_PPM,nb_procedures=nb_proc,type=3), # probabilities of undergoing procedures (columns) in each ward (row), merged for both admission routes; Any time-step is ok but need to be the same as the transition matrix
                 dist_risk=dist_risk, # distribution parameters of the risk for each procedure (dist can be normal, log-normal or pert)
                 nb_material_new=Q_mat_new_init_Low , # quantity of equipment for each type in  each ward
                 nb_material_used=Q_mat_used_init_Low, # quantity of equipment for each type in  each ward
                 nb_material_contaminated=Q_mat_cont_init_Low, # quantity of equipment for each type in  each ward
                 min_e_phase=min_e_phase_HBV, # min value for the eclipse phase
                 max_e_phase=max_e_phase_HBV, # max value for the eclipse phase
                 time_step=time_step, # must be expressed in hours (ex: if time step is 5min then inform 1/12)
                 patient_followup=F, # useless
                 table_proc_equip=as.matrix(corresp_eq_proc), # correspondence table between procedures (col 1) and equipments (col 2)
                 procedure_names=colnames(matrix_proc_eq[,4:18]), # vector of procedure names 
                 equipment_names=matrix_proc_eq$item, # vector of equipments names
                 sterilization_prob=good_ster_Low, # probability of good sterilization for each equipment
                 eq_quantity = "fixed",
                 refill_quantities=refill_quantities_Low, # quantity of equipment to add to the pool of unused equipment for each type at each refill event
                 refill_freq=refill_freq_Low,# refill frequency for each equipment type (must be expressed in hours)
                 equip_bin=round(matrix(nrow = 10,ncol = n_wards,data = runif(280)*100000)),
                 output = "simple",
                 prob_screening=0.72,
                 intervention = intervention,
                 at_risk_wards = at_risk_wards,
                 threshold_reuse= rep(-1,10),
                 pathogen = "HBV",
                 id_sim = n)
}

####################
### INTERVENTIONS ##
####################

################################
### Ward-based interventions ###
################################

prev_init = prev_init_ward
prev_type = "ward"
intervention = "ward-based"

#######################################
### High-resource setting case case ###
#######################################

at_risk_wards = c(28,22,13)

output_WB_High_inter = foreach(n = 1:n_sim ,.packages = c("Rcpp","BloodPaTH"))%dopar%{
  
  set.seed(n)
  
  BloodPaTH_model_inter(t=time, #simulation time
                 n_patients=n_patients, # number of patients
                 nb_wards=28, # number of wards
                 nb_adm=2, # number of admission routes
                 prev_init = prev_init, # initial prevalence for each admission route
                 prev_type = prev_type,
                 WT_matrix=entry_conversion(list_TM,nb_wards=n_wards,type=1), #transition matrix, specific format (merged matrix for all admission routes) but any time-step will do
                 init_prob = rbind(freq.adm.SH,freq.adm.IM), # matrix of probabilities of admission for first ward for both adm routes, nb_adm rows . 
                 adm_prob=c(0.5,0.5),# probabilities of admission for each route
                 nb_procedures=nb_proc, #number of procedures 
                 nb_equipments=10, # number of equipments
                 PPM_matrix=entry_conversion(list_PPM,nb_procedures=nb_proc,type=3), # probabilities of undergoing procedures (columns) in each ward (row), merged for both admission routes; Any time-step is ok but need to be the same as the transition matrix
                 dist_risk=dist_risk, # distribution parameters of the risk for each procedure (dist can be normal, log-normal or pert)
                 nb_material_new=Q_mat_new_init_High , # quantity of equipment for each type in  each ward
                 nb_material_used=Q_mat_used_init_High, # quantity of equipment for each type in  each ward
                 nb_material_contaminated=Q_mat_cont_init_High, # quantity of equipment for each type in  each ward
                 min_e_phase=min_e_phase_HCV, # min value for the eclipse phase
                 max_e_phase=max_e_phase_HCV, # max value for the eclipse phase
                 time_step=time_step, # must be expressed in hours (ex: if time step is 5min then inform 1/12)
                 patient_followup=F, # useless
                 table_proc_equip=as.matrix(corresp_eq_proc), # correspondence table between procedures (col 1) and equipments (col 2)
                 procedure_names=colnames(matrix_proc_eq[,4:18]), # vector of procedure names 
                 equipment_names=matrix_proc_eq$item, # vector of equipments names
                 sterilization_prob=good_ster_High, # probability of good sterilization for each equipment
                 eq_quantity = "fixed",
                 refill_quantities=refill_quantities_High, # quantity of equipment to add to the pool of unused equipment for each type at each refill event
                 refill_freq=refill_freq_High,# refill frequency for each equipment type (must be expressed in hours)
                 equip_bin=matrix(data = 0),
                 output = "simple",
                 prob_screening=1,
                 intervention = intervention,
                 at_risk_wards = at_risk_wards,
                 threshold_reuse= rep(-1,10),
                 id_sim = n)
}

#################################
### Low-resource setting case ###
#################################

at_risk_wards = c(28,27,22)

output_WB_Low_inter = foreach(n = 1:n_sim ,.packages = c("Rcpp","BloodPaTH"))%dopar%{
  
  set.seed(n)
  
  BloodPaTH_model_inter(t=time, #simulation time
                 n_patients=n_patients, # number of patients
                 nb_wards=28, # number of wards
                 nb_adm=2, # number of admission routes
                 prev_init = prev_init, # initial prevalence for each admission route
                 prev_type = prev_type,
                 WT_matrix=entry_conversion(list_TM,nb_wards=n_wards,type=1), #transition matrix, specific format (merged matrix for all admission routes) but any time-step will do
                 init_prob = rbind(freq.adm.SH,freq.adm.IM), # matrix of probabilities of admission for first ward for both adm routes, nb_adm rows . 
                 adm_prob=c(0.5,0.5),# probabilities of admission for each route
                 nb_procedures=nb_proc, #number of procedures 
                 nb_equipments=10, # number of equipments
                 PPM_matrix=entry_conversion(list_PPM,nb_procedures=nb_proc,type=3), # probabilities of undergoing procedures (columns) in each ward (row), merged for both admission routes; Any time-step is ok but need to be the same as the transition matrix
                 dist_risk=dist_risk, # distribution parameters of the risk for each procedure (dist can be normal, log-normal or pert)
                 nb_material_new=Q_mat_new_init_Low , # quantity of equipment for each type in  each ward
                 nb_material_used=Q_mat_used_init_Low, # quantity of equipment for each type in  each ward
                 nb_material_contaminated=Q_mat_cont_init_Low, # quantity of equipment for each type in  each ward
                 min_e_phase=min_e_phase_HCV, # min value for the eclipse phase
                 max_e_phase=max_e_phase_HCV, # max value for the eclipse phase
                 time_step=time_step, # must be expressed in hours (ex: if time step is 5min then inform 1/12)
                 patient_followup=F, # useless
                 table_proc_equip=as.matrix(corresp_eq_proc), # correspondence table between procedures (col 1) and equipments (col 2)
                 procedure_names=colnames(matrix_proc_eq[,4:18]), # vector of procedure names 
                 equipment_names=matrix_proc_eq$item, # vector of equipments names
                 sterilization_prob=good_ster_Low, # probability of good sterilization for each equipment
                 eq_quantity = "fixed",
                 refill_quantities=refill_quantities_Low, # quantity of equipment to add to the pool of unused equipment for each type at each refill event
                 refill_freq=refill_freq_Low,# refill frequency for each equipment type (must be expressed in hours)
                 equip_bin=round(matrix(nrow = 10,ncol = n_wards,data = runif(280)*100000)),
                 output = "simple",
                 prob_screening=1,
                 intervention = intervention,
                 at_risk_wards = at_risk_wards,
                 threshold_reuse= rep(-1,10),
                 pathogen = "HCV",
                 id_sim = n)
}

###################################
### Patient-based interventions ###
###################################

prev_init = prev_init_ward
prev_type = "ward"
intervention = "patient-based"

#######################################
### High-resource setting case case ###
#######################################

at_risk_wards = c(28,27,13)

output_PB_High_inter = foreach(n = 1:n_sim ,.packages = c("Rcpp","BloodPaTH"))%dopar%{
  
  BloodPaTH_model_inter(t=time, #simulation time
                 n_patients=n_patients, # number of patients
                 nb_wards=28, # number of wards
                 nb_adm=2, # number of admission routes
                 prev_init = prev_init, # initial prevalence for each admission route
                 prev_type = prev_type,
                 WT_matrix=entry_conversion(list_TM,nb_wards=n_wards,type=1), #transition matrix, specific format (merged matrix for all admission routes) but any time-step will do
                 init_prob = rbind(freq.adm.SH,freq.adm.IM), # matrix of probabilities of admission for first ward for both adm routes, nb_adm rows . 
                 adm_prob=c(0.5,0.5),# probabilities of admission for each route
                 nb_procedures=nb_proc, #number of procedures 
                 nb_equipments=10, # number of equipments
                 PPM_matrix=entry_conversion(list_PPM,nb_procedures=nb_proc,type=3), # probabilities of undergoing procedures (columns) in each ward (row), merged for both admission routes; Any time-step is ok but need to be the same as the transition matrix
                 dist_risk=dist_risk, # distribution parameters of the risk for each procedure (dist can be normal, log-normal or pert)
                 nb_material_new=Q_mat_new_init_High , # quantity of equipment for each type in  each ward
                 nb_material_used=Q_mat_used_init_High, # quantity of equipment for each type in  each ward
                 nb_material_contaminated=Q_mat_cont_init_High, # quantity of equipment for each type in  each ward
                 min_e_phase=min_e_phase_HCV, # min value for the eclipse phase
                 max_e_phase=max_e_phase_HCV, # max value for the eclipse phase
                 time_step=time_step, # must be expressed in hours (ex: if time step is 5min then inform 1/12)
                 patient_followup=F, # useless
                 table_proc_equip=as.matrix(corresp_eq_proc), # correspondence table between procedures (col 1) and equipments (col 2)
                 procedure_names=colnames(matrix_proc_eq[,4:18]), # vector of procedure names 
                 equipment_names=matrix_proc_eq$item, # vector of equipments names
                 sterilization_prob=good_ster_High, # probability of good sterilization for each equipment
                 eq_quantity = "fixed",
                 refill_quantities=refill_quantities_High, # quantity of equipment to add to the pool of unused equipment for each type at each refill event
                 refill_freq=refill_freq_High,# refill frequency for each equipment type (must be expressed in hours)
                 equip_bin=matrix(data = 0),
                 output = "simple",
                 prob_screening=0.405,
                 intervention = intervention,
                 at_risk_wards = at_risk_wards,
                 threshold_reuse= rep(-1,10),
                 pathogen="HCV",
                 id_sim = n)
}


#################################
### Low-resource setting case ###
#################################

at_risk_wards = c(28,27,13)

output_PB_Low_inter = foreach(n = 1:n_sim ,.packages = c("Rcpp","BloodPaTH"))%dopar%{
  
  BloodPaTH_model_inter(t=time, #simulation time
                 n_patients=n_patients, # number of patients
                 nb_wards=28, # number of wards
                 nb_adm=2, # number of admission routes
                 prev_init = prev_init, # initial prevalence for each admission route
                 prev_type = prev_type,
                 WT_matrix=entry_conversion(list_TM,nb_wards=n_wards,type=1), #transition matrix, specific format (merged matrix for all admission routes) but any time-step will do
                 init_prob = rbind(freq.adm.SH,freq.adm.IM), # matrix of probabilities of admission for first ward for both adm routes, nb_adm rows . 
                 adm_prob=c(0.5,0.5),# probabilities of admission for each route
                 nb_procedures=nb_proc, #number of procedures 
                 nb_equipments=10, # number of equipments
                 PPM_matrix=entry_conversion(list_PPM,nb_procedures=nb_proc,type=3), # probabilities of undergoing procedures (columns) in each ward (row), merged for both admission routes; Any time-step is ok but need to be the same as the transition matrix
                 dist_risk=dist_risk, # distribution parameters of the risk for each procedure (dist can be normal, log-normal or pert)
                 nb_material_new=Q_mat_new_init_Low , # quantity of equipment for each type in  each ward
                 nb_material_used=Q_mat_used_init_Low, # quantity of equipment for each type in  each ward
                 nb_material_contaminated=Q_mat_cont_init_Low, # quantity of equipment for each type in  each ward
                 min_e_phase=min_e_phase_HCV, # min value for the eclipse phase
                 max_e_phase=max_e_phase_HCV, # max value for the eclipse phase
                 time_step=time_step, # must be expressed in hours (ex: if time step is 5min then inform 1/12)
                 patient_followup=F, # useless
                 table_proc_equip=as.matrix(corresp_eq_proc), # correspondence table between procedures (col 1) and equipments (col 2)
                 procedure_names=colnames(matrix_proc_eq[,4:18]), # vector of procedure names 
                 equipment_names=matrix_proc_eq$item, # vector of equipments names
                 sterilization_prob=good_ster_Low, # probability of good sterilization for each equipment
                 eq_quantity = "fixed",
                 refill_quantities=refill_quantities_Low, # quantity of equipment to add to the pool of unused equipment for each type at each refill event
                 refill_freq=refill_freq_Low,# refill frequency for each equipment type (must be expressed in hours)
                 equip_bin=round(matrix(nrow = 10,ncol = n_wards,data = runif(280)*100000)),
                 output = "simple",
                 prob_screening=0.72,
                 intervention = intervention,
                 at_risk_wards = at_risk_wards,
                 threshold_reuse= rep(-1,10),
                 pathogen = "HCV",
                 id_sim = n)
}



# Figure 2 : Daily incidence and ward incidence

a = plot_incidence(output_baseline_High ,50,365*24*12,ylim = c(0,0.4))[[1]]

cuminc_baseline = yearly_cum_inc(output_baseline_High ,50)    

a =  annotate_figure(a, fig.lab = paste0("Yearly cum. inc. : ",cuminc_baseline[1],"/100,000 [",
                                         cuminc_baseline[2],"-",
                                         cuminc_baseline[3],"]"),
                     fig.lab.pos = "top.right",fig.lab.size = 14)


b = plot_inc_wards(output_baseline_High ,50,28,rm_0 = F)[[1]]


c = plot_incidence(output_baseline_Low,50,365*24*12,ylim = c(0,0.4))[[1]]

cuminc_private = yearly_cum_inc(output_baseline_Low  ,50)    

c =  annotate_figure(c, fig.lab = paste0("Yearly cum. inc. : ",cuminc_private[1],"/100,000 [",
                                         cuminc_private[2],"-",
                                         cuminc_private[3],"]"),
                     fig.lab.pos = "top.right",fig.lab.size = 14)

d = plot_inc_wards(output_baseline_Low,50,28,rm_0 = F)[[1]]


plot=ggarrange(a,b,c,d,
               labels = c("A", "B", "C","D"),
               ncol = 2, nrow = 2)


# Figure 3 : Attributable portion for each device 

a = plot_att_eq(baseline_univ_HBV,50,10,plot_log = F)

b = plot_att_eq(baseline_private_HBV,50,10,plot_log = F)

ggarrange(a,b,
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

#Figure 4: Comparing interventions

#Function that compares baseline vs. intervention
comp_inter = function(output,n_sim) {
  
  dif = c()
  for(i in 1:n_sim){
    
    dif[i] = 1- ((sum(output[[i]]$incidence_inter+10^-10))/(sum(output[[i]]$incidence+10^-10))) 
    if(dif[i] < -1) {dif[i]=-1}
  }
  
  
  return(dif)
}

a_inter = comp_inter(output_WB_High_inter,n_sim)

b_inter = comp_inter(output_PB_High_inter,n_sim)

c_inter = comp_inter(output_WB_Low_inter,n_sim)

d_inter = comp_inter(output_PB_Low_inter,n_sim)

#data frame for plotting 
df_plot_inter = data.frame(value = c(a_inter,b_inter,c_inter,d_inter),
                           inter = c(rep("Ward-level",n_sim),rep("Upon-admission",n_sim),rep("Ward-level",n_sim),rep("Upon-admission",n_sim)),
                           Sce = c(rep("High-resource setting",n_sim*2),rep("Low-resource setting",n_sim*2)))

f <- function(x) {
  r <- quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

ggplot(data=df_plot_inter,aes(x=Sce, y=value, fill=inter)) +
  stat_summary(fun.data = f, geom="boxplot", position="dodge")+
  geom_point(aes(fill=inter),position=position_jitterdodge(),alpha=0.3,shape = 21)+
  scale_y_continuous(labels = function(x) paste0(x*100, "%"),limits = c(-1,1))+
  scale_fill_manual(values=c("#009E73","#E69F00"))+
  scale_color_manual(values=c("#009E73","#E69F00"))+
  ylab("Reduction in yearly cumulative incidence (%)")+
  xlab("Setting")+
  theme_bw(base_size = 16)

