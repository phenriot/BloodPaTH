##########################
### Loading packages   ###
##########################
library(BloodPaTH)
library(parallel)
library(graphics)
library(foreach)
library(doParallel)
library(DescTools)
library(reshape2)
library(EnvStats)
library(mc2d)
library(truncnorm)
library(GenSA)
library(ggplot2)
library(ggridges)
library(ggstance)
library(truncdist)
library(dplyr)
library(ggnewscale)
library(paletteer)
library(dichromat)
library(ggthemes)
library(ggtext)
library(lubridate)
library(patchwork)
library(scales)
library(viridis)
library(lvplot)

##############################
### Setting up the cluster ###
##############################

n_cores = detectCores()- 4 # Number of cores allocated to the local cluster

clust <- makeCluster(n_cores)

registerDoParallel(clust)

n_sim = 10

###############
## Functions ##
###############

# Returns a vector of size N summing to 1
rand_vect <- function(N, pos.only = TRUE) {
  
  sd = 1
  prob = rep(0.5,N)
  
  vec <- rbinom(N,1,prob)
  deviation <- 100 - sum(vec)
  for (. in seq_len(abs(deviation))) {
    vec[i] <- vec[i <- sample(N, 1,replace=TRUE,prob)] + sign(deviation)
  }
  if (pos.only) while (any(vec < 0)) {
    negs <- vec < 0
    pos  <- vec > 0
    vec[negs][i] <- vec[negs][i <- sample(sum(negs), 1)] + 1
    vec[pos][i]  <- vec[pos ][i <- sample(sum(pos ), 1)] - 1
  }
  
  return(vec/100)
}

# Returns combined matrices for the lists of transition matrices (type=1) or probability matrices (type=2)
list_to_combined_matrices = function(input,type,nb_wards=NULL,nb_procedures=NULL) {
  
  # For transition matrices 
  if(type==1) {
    
    if(is.null(nb_wards)) {stop("Please provide a value for the nb_wards argument")}
    
    storage_TM = matrix( ncol = nb_wards+1,nrow = 0);
    
    for(i in 1:length(input)) {
      
      storage_TM=rbind(storage_TM,input[[i]])  
      
    } 
    
    return(as.matrix(storage_TM))
    
  }
  
  
  if(type==2) {
    
    if(is.null(nb_procedures)) {stop("Please provide a value for the nb_procedures argument")}
    
    storage_procedures = matrix( ncol = nb_procedures,nrow = 0);
    
    for(i in 1:length(input)) {
      
      storage_procedures =rbind(storage_procedures,input[[i]])  
      
    } 
    
    return(as.matrix(storage_procedures))
  }
  
}

################
## Parameters ##
################


time_step = 1/12 # Expressed in hours, here this corresponds to 5 minutes time steps
time = 8760 # Simulation time, 105,120 5-minutes-time-steps
nb_patients = 970
nb_proc = 15 + 1 # 15 procedures + 1 procedure corresponding to the "no procedure" event 
procedure_names = paste("Procedure",LETTERS[seq( from = 1, to = nb_proc-1 )]) # Vector with the names of each procedure (excluding the "no procedure" event)
nb_devices  = 10 # Number of devices
device_names = paste("Device",LETTERS[seq( from = 1, to = nb_devices )]) # Vector with the names of each device
nb_wards = 28 # Number of wards in the setting
at_risk_wards = c(0) # Vector of IDs of the most at-risk wards
intervention = "none" # Can be "none", "ward-based" or "patient-based"
prob_screening = 0 # Probability of screening, only works when there is an intervention 
nb_adm_routes = 2 # Number of admission routes
adm_prob = c(0.5,0.5) # Admission probability for each route

sterilization_prob = runif(nb_devices,min=0,max = 1) # Random sterilization probabilities for each device

# Load list of transition matrices (list is of length 2 as there is 2 admission routes)
# These matrices are of size (nb_wards + 1)*(nb_wards + 1), the last row and column are the probabilities associated with the event of getting discharged after leaving a particular ward
list_TM = readRDS("C:/Users/paulh/Desktop/Data example/List_Transition_Matrices.rds")
combined_TM = list_to_combined_matrices(list_TM,nb_wards=nb_wards,type=1)

indices_A1 = which(rowSums(list_TM[[1]])[-(nb_wards+1)]==1) # Wards in which patients go within department A1
indices_A2 = which(rowSums(list_TM[[2]])[-(nb_wards+1)]==1) # Wards in which patients go within department A2

freq_adm_A1 = rep(0,28) #Initializing the vector of probabilities of admission for each departments when patient is admitted in department A1
freq_adm_A2 = rep(0,28) #Initializing the vector of probabilities of admission for each departments when patient is admitted in department A1

freq_adm_A1[indices_A1] = rand_vect(length(indices_A1)) # Random probabilities of admission for wards in which patients entering the hospital are admitted in A1
freq_adm_A2[indices_A2] = rand_vect(length(indices_A2)) # Random probabilities of admission for wards in which patients entering the hospital are admitted in A2
#Probabilities of admission for wards within which no patient is hospitalized in a given department need to be set to 0, otherwise error 

# Load list of probability matrices : each element gives the probability of undergoing a particular procedure within a given ward (list is of length 2 as there is 2 admission routes)
# These matrices are of size (nb_wards)*(nb_proc)
list_PPM = readRDS("C:/Users/paulh/Desktop/Data example/List_Proc_Prob_Matrices.rds")
combined_PPM = list_to_combined_matrices(list_PPM,nb_procedures=nb_proc,type=2)

#Matrix of association between devices and procedures
association_matrix = read.csv2("C:/Users/paulh/Desktop/Data example/association_devices_procedures.csv")
association_matrix = as.matrix(association_matrix)

nb_devices_new = matrix(data = round(runif(nb_wards*nb_devices,min = 0,max=100000)),nrow =nb_devices,ncol = nb_wards) # Random matrix for the number of sterile devices in each ward at initialization
nb_devices_used = nb_devices_new*0.2 # Random matrix for the number of non-sterile devices in each ward at initialization
nb_devices_cont = nb_devices_new*0 # We start with not contaminated device at initialization 

refill_time = round(matrix(data = round(runif(nb_wards*nb_devices,min = 0,max=8760)),nrow =nb_devices,ncol = nb_wards)) # Time at which a ward gets refilled with a given type of device, size : nb_wards*nb_devices
refill_quantities = matrix(data = round(runif(nb_wards*nb_devices,min = 0,max=10000)),nrow =nb_devices,ncol = nb_wards) # Associated refill quantities

# Pathogen-associated parameters

initial_ward_prevalence = runif(nb_wards,min = 0, max = 0.05) # Randomly draws prevalence levels for each ward 
prevalence_type = "ward" # Should we consider the prevalence at a setting level of a ward level ? 

# Load a data.frame with the risk of getting infected for each of the procedures
# The first 3 columns correspond to the values of the parameters associated with the distribution of the risk for each of the procedures (rows)
# The last column is the name of the distribution, for now only "lnorm", "norm" and "pert" are working
# The values of the parameters of last row always need to be set to 0 as this corresponds to the "no_procedure" event
dist_risk = read.csv2("C:/Users/paulh/Desktop/Data example/risk_dist.csv")[,3:6]
dist_risk = rbind(dist_risk,c(0,0,NA,"lnorm")) #Adding a "fake" row for the "no procedure" event 
dist_risk$par_1 = as.numeric(dist_risk$par_1)
dist_risk$par_2 = as.numeric(dist_risk$par_2)
dist_risk$par_3 = as.numeric(dist_risk$par_3)


min_e_phase_HCV = 48 # Min. duration of the eclipse phase
max_e_phase_HCV = 336 # Max. duration of the eclipse phase

# Launching the model
output = foreach(n = 1:n_sim ,.packages = c("Rcpp","BloodPaTH"))%dopar%{
  
  BloodPaTH_model(time_step=time_step, # must be expressed in hours (ex: if time step is 5min then inform 1/12)
                  t=time, #simulation time
                  pathogen = "HCV",
                  min_e_phase=min_e_phase_HCV, # min value for the eclipse phase
                  max_e_phase=max_e_phase_HCV, # max value for the eclipse phase
                  nb_patients=nb_patients, # number of patients
                  nb_wards=nb_wards, # number of wards
                  nb_adm=nb_adm_routes, # number of admission routes
                  adm_prob=adm_prob,# probabilities of admission for each route
                  init_prob = rbind(freq_adm_A1,freq_adm_A2), # matrix of probabilities of admission for first ward for both adm routes, nb_adm rows . 
                  prev_type = prevalence_type,
                  prev_init = initial_ward_prevalence, # initial prevalence for each admission route
                  WT_matrix=combined_TM, #transition matrix, specific format (merged matrix for all admission routes) but any time-step will do
                  nb_procedures=nb_proc, #number of procedures 
                  procedure_names=procedure_names, # vector of procedure names 
                  nb_devices=nb_devices, # number of equipments
                  device_names=device_names, # vector of equipments names
                  nb_devices_new=nb_devices_new , # quantity of equipment for each type in  each ward
                  nb_devices_used=nb_devices_used, # quantity of equipment for each type in  each ward
                  nb_devices_contaminated=nb_devices_cont, # quantity of equipment for each type in  each ward
                  refill_quantities=refill_quantities, # quantity of equipment to add to the pool of unused equipment for each type at each refill event
                  refill_freq=refill_time,# refill frequency for each equipment type (must be expressed in hours)
                  table_procedures_devices=association_matrix, # correspondence table between procedures (col 1) and equipments (col 2)
                  sterilization_prob=sterilization_prob, # probability of good sterilization for each equipment        
                  PPM_matrix=combined_PPM, # probabilities of undergoing procedures (columns) in each ward (row), merged for both admission routes; Any time-step is ok but need to be the same as the transition matrix
                  dist_risk=dist_risk, # distribution parameters of the risk for each procedure (dist can be normal, log-normal or pert)
                  intervention = intervention, 
                  prob_screening=prob_screening,
                  at_risk_wards = at_risk_wards,
                  output = "simple",
                  id_sim = n)
}


#Yearly cumulative incidence per 100,000 patients 
yearly_cum_inc(output,n_sim=1)

#Daily incidence over a year
plot_incidence(output,n_sim=1,time)
