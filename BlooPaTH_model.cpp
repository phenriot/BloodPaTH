//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

//[[Rcpp::depends("RcppArmadillo")]]
using namespace Rcpp;
using namespace RcppArmadillo;
using namespace std;

// Loading all needed packages from the R environment
Environment pkg_truncnorm = Environment::namespace_env("truncnorm");
Environment pkg_mc2d = Environment::namespace_env("mc2d");
Environment pkg_stats = Environment::namespace_env("stats");

// Loading all needed function from these packages
Function random_truncnorm=pkg_truncnorm["rtruncnorm"];
Function random_pert = pkg_mc2d["rpert"];
Function random_lnorm= pkg_stats["rlnorm"];
Function random_norm= pkg_stats["rnorm"];
Function random_unif= pkg_stats["runif"];

//Loading a set_seed function to be used within the model
// [[Rcpp::export]]
void set_seed(unsigned int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);  
}

//Data transformation : Returns a list including all associations between procedures and devices to be used within the model's function  
// [[Rcpp::export]]
List transform_list(int nb_proc, arma::imat table_proc_equip) {
  
  List list_proc(nb_proc);
  
  for(int i = 0;i<list_proc.length();i++) {
    
    list_proc(i) =  arma::irowvec(1);
    
  }
  
  for(int i = 0;i<size(table_proc_equip)(0);i++) {
    
    arma::irowvec current_p ;
    current_p = table_proc_equip(i,0);
    arma::irowvec current_e ;
    current_e =table_proc_equip(i,1);
    
    
    for(int j=0; j<list_proc.length();j++ ) {
      
      arma::irowvec storage_vec = list_proc(j);
      
      if(j==current_p(0)-1) {list_proc(j) = join_rows(storage_vec,current_e);}
      
    }
    
  }
  
  for(int i = 0;i<list_proc.length();i++) {
    
    arma::irowvec subset = list_proc(i);
    
    subset.shed_col(0);
    list_proc(i)=subset;
    
  }

  return list_proc;
  
}

///////////////////
// MAIN FUNCTION //
//////////////////

// Returns a list of elements : 
// incidence : a vector of length t giving the number of new cases at each time-step at a setting level
// new_patients : a vector of length t giving the number of new patients entering the setting at each time-step
// s_patients : a vector of length t giving the number of susceptible patients in the setting at each time-step
// cont_new_patients : a vector of length t giving the number of newly contaminated patients in the setting at each time-step
// ward_s : a matrix of size nb_wards*t giving the number of susceptible patients at a ward level at each time-step
// ward_event : a vector of length nb_wards giving the total number of contamination events in the setting across the period
// count_patients : an integer reporting the total number of patients that entered the setting across the period 
// cont_mat : a vector of length nb_devices giving the total number of contaminated devices at the end of the period
// used_mat : a vector of length nb_devices giving the total number of non-sterile devices at the end of the period
// list_proc_equip : a list of length nb_procedures; each element gives the type of devices (coded as digits) used in each procedure
// vec_nb_proc : a vector of length nb_procedures reporting the total number of procedures performed on patients within the setting across the period
// inf_eq : a vector of length nb_procedures reporting the total number of contamination events that occurred after each type of procedure across the period
// count_patient_ward : a vector  providing the number of accesses (not unique accesses) in each ward across the period
// count_tests : an integer reporting the total number of tests performed across the period 
// eq_usage : a matrix of size nb_devices*nb_wards reporting the total effective usage of each device type in each ward across the period


// [[Rcpp::export]]
List BlooPaTH_model(int t, //Simulation time
                    int n_patients, //Number of patients
                    int nb_wards, //Number of wards
                    int nb_adm, //Number of admission routes
                    arma::drowvec prev_init, //Initial prevalence for each admission route
                    std::string prev_type, //"Ward" level or "hospital" level
                    arma::dmat WT_matrix, //Transition matrix, specific format (merged matrix for all admission routes) but any time-step will do
                    arma::dmat init_prob, //Matrix of probabilities of admission for first ward for both admission routes, nb_adm rows .
                    arma::dcolvec adm_prob,//Probabilities of admission for each route
                    int nb_procedures, //Number of procedures (must include the "no procedure" procedure)
                    int nb_devices, //Number of devices
                    arma::dmat PPM_matrix, //Probabilities of undergoing procedures (columns) in each ward (row), merged for both admission routes; Any time-step is ok but need to be the same as the transition matrix
                    Rcpp::DataFrame dist_risk, //Distribution parameters of the risk for each procedure (dist can be normal, log-normal or pert)
                    arma::imat nb_material_new, //Quantity of equipment for each type in  each ward
                    arma::imat nb_material_used, //Quantity of equipment for each type in  each ward
                    arma::imat nb_material_contaminated, //Quantity of equipment for each type in  each ward
                    int min_e_phase, // Min. value for the eclipse phase
                    int max_e_phase, // Max. value for the eclipse phase
                    float time_step, // Must be expressed in hours (ex: if time step is 5min then inform 1/12)
                    bool patient_followup, // NOT USED IN THIS VERSION
                    arma::imat table_proc_equip, // Correspondence table between procedures (column 1) and devices (column 2)
                    StringVector procedure_names, // Vector of procedures names
                    StringVector equipment_names, // Vector of devices names
                    arma::drowvec sterilization_prob, //Vector of probabilities of efficient sterilization for each device
                    std::string eq_quantity, //Should the quantity be fixed ? "fixed" (then refill quantities not used) or "variable" - SECOND OPTION NOT DEFINED IN THIS VERSION
                    arma::imat refill_quantities, // Quantity of devices to add to the pool of unused equipment for each type at each refill event
                    arma::fmat refill_freq,// Refill frequency for each equipment type in each ward (must be expressed in hours)
                    arma::fmat equip_bin,// NOT USED IN THIS VERSION
                    std::string output="simple",// Only one type of output for now
                    float prob_screening=0, // Probability of systematic screening upon admission
                    std::string intervention = "none",//none, patient-based or ward-based
                    arma::irowvec at_risk_wards=0, // vector giving the most at risk wards
                    arma::frowvec threshold_reuse=0, // Threshold before equipment reuse: NOT USED IN THIS VERSION
                    std::string pathogen = "HCV", // Pathogen : needs to be "HCV" or "HBV"
                    int id_sim=1) //Simulation ID
{
  
  List list_proc_equip = transform_list(nb_procedures,table_proc_equip);
  
  
  //clean freq en vecteur pour associer a chaque type de materiel
  arma::dmat storage_adm_TM;
  arma::dmat storage_adm_PPM;
  arma::drowvec storage_init_prob;
  
  StringVector subset_col_proc = table_proc_equip[0];
  StringVector subset_col_equip = table_proc_equip[0];
  
  //vector of indices to recover each individual transition matrix as the entered matrix is of size ((nb_wards+1)*nb_adm) x (nb_wards+1)
  arma::irowvec index_TM = arma::linspace<arma::irowvec>(0,(nb_wards+1)*nb_adm,nb_adm+1);
  
  arma::irowvec index_PPM = arma::linspace<arma::irowvec>(0,(nb_wards)*nb_adm,nb_adm+1);
  
  // error messages
  if(prev_type=="ward" && prev_init.size() != nb_wards) {stop("Error: Initial prevalence type chosen is at a ward-level, should be the same size as the number of wards");}
  if(prev_type=="hospital" && prev_init.size() != nb_adm) {stop("Error: Initial prevalence type chosen is at a hospital-level, should be the same size as the number of admission routes");}
  
  //Vectors used within the function for the epidemiological process
  arma::drowvec pos_neg_vec = {0,2}; // vector with 2 elements : 0 = negative - 2 = positive 
  arma::dcolvec cont_neg_vec = {0,1}; // vector with 2 elements : 0 = negative - 1 = exposed
  arma::dcolvec failure_success = {0,1}; // vector with 2 elements : 0 = failure - 1 = success
  arma::drowvec screened_pos_neg_vec = {1,2}; // vector with 2 elements : 0 = negative - 1 = exposed
  
  
  //number of patients in the hospital
  int number_pat = n_patients;
  
  //Matrix for status of patients in the hospital
  arma::dmat pop_hosp_status = arma::dmat(n_patients,t+1,arma::fill::zeros); 
  
  // counting nb of patients hospitalized on the periode 
  int count_patients = 0;
  
  // counting nb of tests
  int count_tests = 0;
  
  // counting nb of patients hospitalized on the period in each ward
  arma::irowvec count_patient_ward(nb_wards,arma::fill::zeros);
  
  //Matrix for location of patients in the hospital
  arma::dmat pop_hosp_loc = arma::dmat(n_patients,t+1,arma::fill::zeros); 
  
  //Matrix of procedures undergone by patients during hospitalization
  arma::dmat pop_proc = arma::dmat(n_patients,t+1,arma::fill::zeros); 
  
  //Matrix of number of previously used materials
  arma::imat Npu_mat = nb_material_used;
  
  //Matrix of number of available material
  arma::imat N_mat = nb_material_new;
  
  //Matrix of number of contaminated materials for each ward and each procedure (C compartment)
  arma::imat C_mat = nb_material_contaminated;
  
  //Matrix storing all departures of patients
  arma::umat patient_departure;
  
  //Sequence for wards : 1 --> nb_wards
  arma::drowvec wards_seq = arma::linspace<arma::drowvec>(1,nb_wards,nb_wards);
  
  //Sequence for admission places 1 --> nb_adm
  arma::drowvec adm_seq = arma::linspace<arma::drowvec>(1,nb_adm,nb_adm);
  
  //Sequence for admission places 1 --> nb_procedures
  arma::drowvec proc_seq = arma::linspace<arma::drowvec>(1,nb_procedures,nb_procedures);
  
  //Incidence vector
  arma::drowvec incidence(t,arma::fill::zeros);
  
  //Number of susceptible patients vector
  arma::drowvec s_patients(t,arma::fill::zeros);
  
  //New patients within the setting
  arma::drowvec new_patient_vec(t,arma::fill::zeros);
  new_patient_vec(0) = n_patients; // initialization at t0
  
  //Contaminated new patients
  arma::drowvec cont_new_patient_vec(t,arma::fill::zeros);
  
  //Time counter to assign patient in the E state to the I state
  arma::drowvec state_counter(n_patients,arma::fill::zeros);// starts at 0
  
  //Refill counter for devices
  arma::imat refill_counter(nb_devices,nb_wards,arma::fill::zeros);// starts at 0
  
  //Time counter to tell when used devices need to be thrown away (not used in baseline scenarios)
  arma::imat bin_counter(nb_devices,nb_wards,arma::fill::zeros);// starts at 0
  
  //Matrix of devices usage
  arma::imat eq_usage(nb_devices,nb_wards,arma::fill::zeros);// starts at 0
  
  //Vector storing admission route for each patient
  arma::drowvec admission_route(n_patients,arma::fill::zeros);// starts at 0
  
  //Vector storing screening status for each patient 
  arma::drowvec screened(n_patients,arma::fill::zeros);// starts at 0
  
  //Vector storing the probability of getting screened
  arma::dcolvec screen_prob_vec = {1-prob_screening,prob_screening};
  
  //Matrix to know in patient has been screened for HCV or HBV or HIV , 0 = not screened, 1 = screened and - , 2 = screened and +
  arma::imat ward_screening = arma::imat(n_patients,nb_wards,arma::fill::zeros); //No need to t+1 because first index = 0
  
  //Vector to store nb of infection events in each ward
  arma::drowvec ward_event(nb_wards,arma::fill::zeros);// starts at 0
  
  //Matrix to store nb of susceptible patients in each ward at each time step
  arma::dmat ward_s = arma::dmat(nb_wards,t+1,arma::fill::zeros);
  
  //vector to store nb of infection events
  arma::drowvec cont_mat(nb_devices,arma::fill::zeros);// starts at 0
  
  //vector to store nb of infection events
  arma::drowvec used_mat(nb_devices,arma::fill::zeros);// starts at 0
  
  //vector to store nb of proc by proc
  arma::drowvec vec_nb_proc(nb_procedures,arma::fill::zeros);// starts at 0
  
  //vector to store nb of infection per equipment type
  arma::drowvec inf_eq(nb_devices,arma::fill::zeros);// starts at 0
  
  //Matrix to know if a patient has already been into a given ward
  arma::imat mat_patient_ward = arma::imat(n_patients,nb_wards,arma::fill::zeros); //No need to t+1 because first index = 0
  
  //A vector of zeros used when a new patient enters the setting : initializes to zero a patient's ward history
  arma::irowvec init_patient_wards(nb_wards,arma::fill::zeros);
  
  
  //Initializes the eclipse phase value for each patient
  arma::drowvec e_phase = Rcpp::runif( n_patients, min_e_phase, max_e_phase);
  e_phase = round(e_phase);
  
  //Vector storing all patients IDs
  arma::drowvec patient_id =  arma::linspace<arma::drowvec>(0,n_patients-1,n_patients);
  
  //Loop simulating the initial patients status and screening policies upon admission for all patients (at t0) 
  for(int i = 0; i<n_patients;i++) {
    
    admission_route(i) = sample(adm_seq,1,false,adm_prob)(0);
    
    storage_init_prob= init_prob.row(admission_route(i)-1);
    
    pop_hosp_loc(i,0) = sample(wards_seq,1,false,storage_init_prob.t())(0);
    
    if(prev_type == "hospital") {
      
      arma::dcolvec pos_neg_prob = {1-prev_init((admission_route(i)-1)),prev_init((admission_route(i)-1))};
      
      // Is patient i infected ? Depends on admission route
      pop_hosp_status(i,0) = sample(pos_neg_vec,1,false,pos_neg_prob)(0);
      
      if(pop_hosp_status(i,0)==2) {cont_new_patient_vec(0) = cont_new_patient_vec(0)+1;}
      
    }
    
    if(prev_type == "ward") {
      
      arma::dcolvec pos_neg_prob = {1-prev_init((pop_hosp_loc(i,0)-1)),prev_init((pop_hosp_loc(i,0)-1))};
      
      // Is patient i infected ? Depends on admission route
      pop_hosp_status(i,0) = sample(pos_neg_vec,1,false,pos_neg_prob)(0);
      
      if(pop_hosp_status(i,0)==2) {cont_new_patient_vec(0) = cont_new_patient_vec(0)+1;}
      
    }
    
    
    if(intervention == "ward-based" && std::find(at_risk_wards.begin(), at_risk_wards.end(), (pop_hosp_loc(i,0))) != at_risk_wards.end()) {
      
      int screen = sample(failure_success,1,false,screen_prob_vec)(0);
      
      if (screen==1 && pop_hosp_status(i,0)==0){screened(i)=1;  count_tests = count_tests +1 ;}
      if (screen==1 && pop_hosp_status(i,0)==2){screened(i)=2;  count_tests = count_tests +1 ;}
      
    }
    
    if(intervention =="patient-based") {
      
      int screen = sample(failure_success,1,false,screen_prob_vec)(0);
      
      if (screen==1 && pop_hosp_status(i,0)==0) {screened(i)=1; count_tests = count_tests +1 ;}
      if (screen==1 && pop_hosp_status(i,0)==2) {screened(i)=2; count_tests = count_tests +1 ;}
      
    }
  }
  
  
  
  //Adding a supplementary state for patients leaving the hospital (state "out")
  int state_out = nb_wards+1;
  arma::drowvec wards_out_seq = wards_seq;
  wards_out_seq.resize(state_out);
  wards_out_seq(wards_out_seq.size()-1) = state_out;
  
  //Loop simulating the main process  
  for(unsigned int time = 0; time<t;time++) {
    
    rdm_patient  = sample(patient_id,n_patients,false); //randomly sorting patients ID
    
    // Looping over all patients 
    for(unsigned int j = 0; j<n_patients;j++) {
      
      unsigned int p = rdm_patient(j); // ID of patient j 
      
      int current_ward = pop_hosp_loc(p,time); // current ward ID
      
      // Storing transition matrix corresponding to a given admission route into an object : rows in span  index_TM(0) --> index_TM(1)-1
      storage_adm_TM = WT_matrix(arma::span(index_TM(admission_route(p)-1),index_TM(admission_route(p))-1),arma::span(0,nb_wards));
      storage_adm_TM =  storage_adm_TM.t();
      arma::dcolvec prob_ward =storage_adm_TM.col(current_ward-1);//minus 1 because first index = 0
      
      //Same thing for the PPM matrix
      storage_adm_PPM = PPM_matrix(arma::span(index_PPM(admission_route(p)-1),index_PPM(admission_route(p))-1),arma::span(0,(nb_procedures-1)));
      storage_adm_PPM =storage_adm_PPM.t();
      
      pop_proc(p,time) = sample(proc_seq,1,false, storage_adm_PPM.col(pop_hosp_loc(p,time)-1))(0); // random draw of proc in prob linked to selected ward
      
      int current_proc = pop_proc(p,time) ; // current procedure the patient undergoes
      
      int viremic_status = pop_hosp_status(p,time) ; // viremic status of patient p
      
      
      //Adding 1 to the passage history for ward which ID=current_ward
      if(mat_patient_ward(p,current_ward-1) == 0) {
        
        count_patient_ward(current_ward-1) = count_patient_ward(current_ward-1)+1;
        
      }
      
      //Stores the fact that patient p has been in ward current_ward
      mat_patient_ward(p,current_ward-1) = 1;
      
      //Adding 1 to the vector counting the number of procedures 
      vec_nb_proc(current_proc-1) =vec_nb_proc(current_proc-1)+1;
      
      if(time>0) {
        
        if((intervention == "ward-based") && 
           screened(p) == 0 && 
           std::find(at_risk_wards.begin(), at_risk_wards.end(), current_ward) != at_risk_wards.end() &&
           pop_hosp_loc(p,time-1) != pop_hosp_loc(p,time)) {
          
          //Random draw for patient's screening
          int screen = sample(failure_success,1,false,screen_prob_vec)(0);
          
          if (screen==1 && pop_hosp_status(p,time-1)==0) {screened(p) = 1;  count_tests = count_tests +1 ;} // patient gets screened and is negative
          if (screen==1 && pop_hosp_status(p,time-1)==2) {screened(p) = 2;  count_tests = count_tests +1 ;} // patient gets screened and is positive
          
        } 
      }
      
      ////////////////////////
      // INFECTED PATIENTS //
      ///////////////////////
      
      // NB :current_proc != nb_proc because the nb_procth procedure refers to "no procedure"
      if((viremic_status == 2 or viremic_status == 3) && current_proc != nb_procedures && current_ward != state_out) {
        
        
        //List of involved devices for a given procedure
        arma::irowvec list_equip = list_proc_equip(current_proc-1);
        
        //Loop over these devices
        for(int eq = 0; eq< list_equip.n_cols; eq++ ) {
          
          int current_device = list_equip(eq); 
          
          eq_usage(current_device-1,current_ward-1) = eq_usage(current_device-1,current_ward-1) + 1;
          
          float n_available = N_mat(current_device-1,current_ward-1); //nb of total available devices for a given procedure and ward
          
          float n_pu = Npu_mat(current_device-1,current_ward-1); //nb of total previously used devices for a given procedure and ward
          
          
          //Initializing an element storing the "reuse" event
          int reuse = 0;
          
          //Total number of available devices (sterile + non-sterile)
          float init_eq = n_available + n_pu;
          
          //If there is no more available device
          if(n_available < 1) {reuse = 1;}
          
          //If reusing device random draw to know if this device is contaminated
          if(reuse == 1) {
            
            used_mat(current_device-1) = used_mat(current_device-1)+1;
            
            float n_c = C_mat(current_device-1,current_ward-1);
            double prob_Cpw = n_c / n_pu; // Probability of drawing a contaminated device
            
            if(prob_Cpw > 1) {prob_Cpw=1;}
            
            arma::dcolvec vec_prob_Cpw = {1-prob_Cpw,prob_Cpw};
            
            int draw_mat = sample(failure_success,1,false,vec_prob_Cpw)(0);
            
            //If device was not previously contaminated, it becomes contaminated with a probability equal to sterilization efficiency for the given device
            if(draw_mat==0){
              
              // Probability of device being well sterilized
              double prob_ster = sterilization_prob(current_device-1);
              arma::dcolvec vec_prob_ster = {1-prob_ster,prob_ster};
              int sterilized = sample(failure_success,1,false,vec_prob_ster)(0);
              
              //If the patient is known as positive then sterilization is 100% efficient (only used when assessing interventions)
              if( (intervention == "patient-based" or intervention == "ward-based") && screened(p) == 2) {sterilized = 1;}
              
              //If sterilization fails
              if(sterilized==0) {C_mat(current_device-1,current_ward-1) = C_mat(current_device-1,current_ward-1)+1;}
              
            }
            
            //If the device was previously contaminated, it can be sterilized
            if(draw_mat==1){
              
              cont_mat(current_device-1) = cont_mat(current_device-1)+1;
              
              // Probability of device being well sterilized
              double prob_ster = sterilization_prob(current_device-1);
              arma::dcolvec vec_prob_ster = {1-prob_ster,prob_ster};
              int sterilized = sample(failure_success,1,false,vec_prob_ster)(0);
              
              //If the patient is known as positive then sterilization is 100% efficient (only used when assessing interventions)
              if( (intervention == "patient-based" or intervention == "ward-based") && screened(p) == 2) {sterilized = 1;}
              
              //If sterilization is effective
              if(sterilized==1) {C_mat(current_device-1,current_ward-1) = C_mat(current_device-1,current_ward-1)-1;}
              
            }
          }
          
          //If the used device is sterile then it gets contaminated depending on the sterilization probability
          if(reuse == 0) {
            
            //Adding 1 to matrix of previously used material
            Npu_mat(current_device-1,current_ward-1) = Npu_mat(current_device-1,current_ward-1) + 1;
            
            //Removing 1 to matrix of available material
            N_mat(current_device-1,current_ward-1) = N_mat(current_device-1,current_ward-1)-1;
            
            // Probability of device being well sterilized
            double prob_ster = sterilization_prob(current_device-1);
            arma::dcolvec vec_prob_ster = {1-prob_ster,prob_ster};
            int sterilized = sample(failure_success,1,false,vec_prob_ster)(0);
            
            //If the patient is known as positive then sterilization is 100% efficient (only used when assessing interventions)
            if( (intervention == "patient-based" or intervention == "ward-based") && screened(p) == 2) {sterilized = 1;}
            
            //If sterilization fails
            if(sterilized==0) {C_mat(current_device-1,current_ward-1) = C_mat(current_device-1,current_ward-1)+1;}
            
          }
        }
      }
      
      //////////////////////
      // EXPOSED PATIENTS //
      //////////////////////
      
      // NB :current_proc != nb_proc because the nb_procth procedure refers to "no procedure"
      if (viremic_status == 1 && current_proc != nb_procedures && current_ward != state_out) {
        
        //List of involved devices for a given procedure
        arma::irowvec list_equip = list_proc_equip(current_proc-1);
        
        //Loop over these devices
        for(int eq = 0; eq< list_equip.n_cols; eq++ ) {
          
          int current_device = list_equip(eq);
          
          eq_usage(current_device-1,current_ward-1) = eq_usage(current_device-1,current_ward-1) + 1;
          
          float n_available = N_mat(current_device-1,current_ward-1); //nb of total available devices for a given proc and ward
          
          float n_pu = Npu_mat(current_device-1,current_ward-1); //nb of total previously used devices for a given procedure and ward
          
          
          //Initializing an element storing the reuse event
          int reuse = 0;
          
          //Total number of available devices (sterile + non-sterile)
          float init_eq = n_available + n_pu;
          
          //If there is no more available device
          if(n_available < 1) {reuse = 1;}
          
          //If there is no reuse
          if(reuse == 0) {
            
            Npu_mat(current_device-1,current_ward-1) = Npu_mat(current_device-1,current_ward-1) + 1;
            N_mat(current_device-1,current_ward-1) = N_mat(current_device-1,current_ward-1)-1;
            
          }
          
          //If reusing device random draw to know if this device is contaminated
          if(reuse == 1) {
            
            used_mat(current_device-1) = used_mat(current_device-1)+1;
            
            float n_c = C_mat(current_device-1,current_ward-1);
            double prob_Cpw = n_c / n_pu; // Probability of drawing a contaminated device
            
            if(prob_Cpw > 1) {Rcout<< prob_Cpw; prob_Cpw=1;}
            
            arma::dcolvec vec_prob_Cpw = {1-prob_Cpw,prob_Cpw};
            
            int draw_mat = sample(failure_success,1,false,vec_prob_Cpw)(0);
            
            //If device was not previously contaminated, it becomes contaminated with a probability equal to sterilization efficiency for the given device
            if(draw_mat==1){
              
              cont_mat(current_device-1) = cont_mat(current_device-1)+1;
              
              // Probability of device being well sterilized
              double prob_ster = sterilization_prob(current_device-1);
              arma::dcolvec vec_prob_ster = {1-prob_ster,prob_ster};
              int sterilized = sample(failure_success,1,false,vec_prob_ster)(0);
              
              //If sterilization is effective
              if(sterilized==1) {C_mat(current_device-1,current_ward-1) = C_mat(current_device-1,current_ward-1)-1;}
              
            }
          }
        }
      }
      
      ///////////////////////////
      // SUSCEPTIBLE PATIENTS //
      //////////////////////////
      
      // NB :current_proc != nb_proc because the nb_procth procedure refers to "no procedure"
      if (viremic_status == 0 && current_proc != nb_procedures && current_ward != state_out) {
        
        //List of involved devices for a given procedure
        arma::irowvec list_equip = list_proc_equip(current_proc-1);
        
        // Vector of the probabilties of devices being contaminated (used if multiple devices for one procedure)
        arma::drowvec prob_eq_vec(list_equip.n_cols);
        
        //Loop over these devices
        for(int eq = 0; eq< list_equip.n_cols; eq++ ) {
          
          int current_device = list_equip(eq); 
          
          eq_usage(current_device-1,current_ward-1) = eq_usage(current_device-1,current_ward-1) + 1;
          
          float n_available = N_mat(current_device-1,current_ward-1); //nb of total available material for a given proc and ward
          
          float n_pu = Npu_mat(current_device-1,current_ward-1); //nb of total previously used devices for a given procedure and ward
          
          //Initializing an element storing the reuse event
          int reuse = 0;
          
          //Total number of available devices (sterile + non-sterile)
          float init_eq = n_available + n_pu;
          
          //If there is no more available device
          if(n_available < 1) {reuse = 1;}
          
          //If reusing device random draw to know if this device is contaminated
          if(reuse == 1) {
            
            used_mat(current_device-1) = used_mat(current_device-1)+1;
            
            float n_c = C_mat(current_device-1,current_ward-1);
            double prob_Cpw = n_c / n_pu; // Probability of drawing a contaminated device
            
            if(prob_Cpw > 1) {Rcout<< prob_Cpw; prob_Cpw=1;}
            
            arma::dcolvec vec_prob_Cpw = {1-prob_Cpw,prob_Cpw};
            
            int draw_mat = sample(failure_success,1,false,vec_prob_Cpw)(0);
            
            prob_eq_vec(eq) = draw_mat;
            
            //If the device was previously contaminated, it can be sterilized
            if(draw_mat==1){
              
              cont_mat(current_device-1) = cont_mat(current_device-1)+1;
              
              // Probability of device being well sterilized
              double prob_ster = sterilization_prob(current_device-1);
              arma::dcolvec vec_prob_ster = {1-prob_ster,prob_ster};
              int sterilized = sample(failure_success,1,false,vec_prob_ster)(0);
              
              //If sterilization is effective
              if(sterilized==1) {C_mat(current_device-1,current_ward-1) = C_mat(current_device-1,current_ward-1)-1;}
              
            }
          }
          
          //If there is no reuse
          if(reuse == 0) {
            
            Npu_mat(current_device-1,current_ward-1) = Npu_mat(current_device-1,current_ward-1) + 1;
            N_mat(current_device-1,current_ward-1) = N_mat(current_device-1,current_ward-1)-1;
            
          }
        }
        
        //If there is at least 1 contaminated device
        if(sum(prob_eq_vec) >= 1) {
          
          StringVector vec_dist = dist_risk["dist"]; // Reaching the type of distribution  
          
          // If the distribution of the risk is log-normal
          if(vec_dist(current_proc-1)=="lnorm") {
            
            arma::fvec lmean_vec = dist_risk[0];
            
            arma::fvec lsd_vec = dist_risk[1];
            
            float lmean_risk = lmean_vec(current_proc-1); //Parameter 1 for the log-normal distribution
            
            float lsd_risk = lsd_vec(current_proc-1); //Parameter 2 for the log-normal distribution
            
            float risk_p = Rcpp::as<float> (random_lnorm(_["n"]=1,_["meanlog"]=lmean_risk,_["sdlog"]=lsd_risk)); //Risk of getting HCV contaminated if the device is contaminated
            
            float prob_E =risk_p ; 
            
            
            // HVB case
            if(pathogen =="HBV") {
              
              float HBV_mult = Rcpp::as<float> (random_unif(_["n"]=1,_["min"]=4,_["max"]=20)); //Risk multiplier for the HBV case
              
              prob_E = risk_p*HBV_mult;
              
              Rcout<< prob_E;
              
            }
            
            //Probability of going from the Susceptible to the Exposed state
            arma::dcolvec prob_inf_vector = {1-prob_E,prob_E};
            
            int draw_infection = sample(cont_neg_vec,1,false,prob_inf_vector)(0);
            
            if (draw_infection == 1) {pop_hosp_status(p,time) = 1; incidence(time)=incidence(time)+1; ward_event(current_ward-1)= ward_event(current_ward-1)+1; 
            
            for(int ei = 0; ei < list_equip.n_cols; ei++) {if(prob_eq_vec(ei)==1) {inf_eq(list_equip(ei)-1) = inf_eq(list_equip(ei)-1)+1;} } 
            
            }
          }
          
          // If the distribution of the risk is gaussian
          if(vec_dist(current_proc-1)=="norm") {
            
            arma::fvec mean_vec = dist_risk[0];
            
            arma::fvec sd_vec = dist_risk[1];
            
            float mean_risk = mean_vec(current_proc-1); //Parameter 1 for the gaussian distribution (mean)
            
            float sd_risk = sd_vec(current_proc-1); //Parameter 2 for the gaussian distribution (sd)
            
            float risk_p = Rcpp::as<float> (random_norm(_["n"]=1,_["mean"]=mean_risk,_["sd"]=sd_risk)); //risk of getting HCV contaminated if material is contaminated
            
            float prob_E = risk_p;
            
            // HVB case
            if(pathogen=="HBV") {
              
              float HBV_mult = Rcpp::as<float> (random_unif(_["n"]=1,_["min"]=4,_["max"]=20)); //risk mult for HBV
              
              prob_E = risk_p*HBV_mult;
              
            }
            
            //Probability of going from the Susceptible to the Exposed state
            arma::dcolvec prob_inf_vector = {1-prob_E,prob_E};
            
            int draw_infection = sample(cont_neg_vec,1,false,prob_inf_vector)(0);
            
            if (draw_infection == 1) {pop_hosp_status(p,time) = 1; incidence(time)=incidence(time)+1; ward_event(current_ward-1)= ward_event(current_ward-1)+1;
            
            for(int ei = 0; ei < list_equip.n_cols; ei++) {if(prob_eq_vec(ei)==1) {inf_eq(list_equip(ei)) = inf_eq(list_equip(ei)-1)+1;} } 
            
            }
          }
          
          // If the risk is pert-distributed
          if(vec_dist(current_proc-1)=="pert") {
            
            arma::fvec min_vec = dist_risk[0];
            
            arma::fvec mode_vec = dist_risk[1];
            
            arma::fvec max_vec = dist_risk[2];
            
            float min_risk = min_vec(current_proc-1); //Parameter 1 for the PERT distribution
            
            float mode_risk = mode_vec(current_proc-1); //Parameter 2 for the PERT distribution
            
            float max_risk = max_vec(current_proc-1); //Parameter 3 for the PERT distribution
            
            float risk_p = Rcpp::as<float> (random_pert(_["n"]=1,_["min"]=min_risk,_["mode"]=mode_risk,_["max"]=max_risk))/100; //risk of getting HCV contaminated if material is contaminated
            
            float prob_E = risk_p;
            
            // HVB case
            if(pathogen =="HBV") {
              
              float HBV_mult = Rcpp::as<float> (random_unif(_["n"]=1,_["min"]=4,_["max"]=20)); //risk mult for HBV
              
              prob_E = risk_p*HBV_mult;
              
            }
            
            //Probability of going from the Susceptible to the Exposed state
            arma::dcolvec prob_inf_vector = {1-prob_E,prob_E};
            
            int draw_infection = sample(cont_neg_vec,1,false,prob_inf_vector)(0);
            
            if (draw_infection == 1) {pop_hosp_status(p,time) = 1; incidence(time)=incidence(time)+1; ward_event(current_ward-1)= ward_event(current_ward-1)+1;
            
            for(int ei = 0; ei < list_equip.n_cols; ei++) {if(prob_eq_vec(ei)==1) {inf_eq(list_equip(ei)-1) = inf_eq(list_equip(ei)-1)+1;} } 
            
            }
          }
        }
      }
      
      //Adding 1 to the status counter if patient is Exposed
      if(pop_hosp_status(p,time) == 1) {state_counter(p) = state_counter(p) + 1;}
      
      //If the counter exceed the duration of the eclipse phase then patient goes from the Exposed to the Infected state
      if(state_counter(p) >= e_phase(p) && pop_hosp_status(p,time) == 1 ) {pop_hosp_status(p,time) = 3;}
      
      pop_hosp_loc(p,time+1) = sample(wards_out_seq,1,false,prob_ward)(0); //Randomly drawing next ward for patient p
      
      pop_hosp_status(p,time+1) = pop_hosp_status(p,time); //Storing patient status at time t
      
      if(pop_hosp_status(p,time)==0) {s_patients(time)=s_patients(time)+1; ward_s(current_ward-1,time)=ward_s(current_ward-1,time)+1;}
      
      // If a patient leaves the hospital, another one enters
      if( pop_hosp_loc(p,time+1) == state_out) {
        
        //Initialize to zero so that history of passing through wards is clear for new patient
        mat_patient_ward.row(p) = init_patient_wards;
        
        //Same thing for screening history
        ward_screening.row(p) = init_patient_wards;
        
        //Add 1 to patients counter
        count_patients = count_patients + 1;
        
        new_patient_vec(time) = new_patient_vec(time) + 1;
        
        
        arma::ucolvec new_departure = {time+1,p};
        patient_departure = join_horiz(patient_departure ,new_departure) ;
        
        admission_route(p) = sample(adm_seq,1,false,adm_prob)(0);
        storage_init_prob= init_prob.row(admission_route(p)-1);
        pop_hosp_loc(p,time+1) = sample(wards_seq,1,false,storage_init_prob.t())(0);//if patient leaves the hospital, another patient instantly replaces him
        
        state_counter(p) = 0; //If patient p leaves the hospital, the counter of infection for patient p goes back to 0
        
        if(prev_type == "ward") {
          
          arma::dcolvec pos_neg_prob = {1-prev_init((pop_hosp_loc(p,time+1)-1)),prev_init((pop_hosp_loc(p,time+1)-1))};
          
          //Is patient p infected ? Depends on admission route
          pop_hosp_status(p,time+1) = sample(pos_neg_vec,1,false,pos_neg_prob)(0);
          
        }
        
        if(prev_type == "hospital") {
          
          arma::dcolvec pos_neg_prob = {1-prev_init((admission_route(p)-1)),prev_init((admission_route(p)-1))};
          
          //Is patient p infected ? Depends on admission route
          pop_hosp_status(p,time+1) = sample(pos_neg_vec,1,false,pos_neg_prob)(0);
          
        }
        
        //Randomly draws an eclipse phase duration
        e_phase(p) = round(Rcpp::runif(1,min_e_phase,max_e_phase)(0));
        
        screened(p)=0;
        
        if(intervention == "ward-based" && std::find(at_risk_wards.begin(), at_risk_wards.end(), (pop_hosp_loc(p,time+1))) != at_risk_wards.end()) {
          
          int screen = sample(failure_success,1,false,screen_prob_vec)(0);
          
          if (screen==1 && pop_hosp_status(p,time+1)==0) {screened(p)=1; count_tests = count_tests +1 ;}
          if (screen==1 && pop_hosp_status(p,time+1)==2) {screened(p)=2; count_tests = count_tests +1 ;}
          
        }
        
        if(intervention =="patient-based") {
          
          int screen_new = sample(failure_success,1,false,screen_prob_vec)(0);
          
          // Rcout << "ici 3 ";
          
          if (screen_new==1 && pop_hosp_status(p,time+1)==0) {screened(p)=1;  count_tests = count_tests +1 ;}
          
          if (screen_new==1 && pop_hosp_status(p,time+1)==2) {screened(p)=2;  count_tests = count_tests +1 ;}
          
        }
      }
    }
    
    //Check if there is a need to refill or throw away devices
    for(unsigned int rf=0;rf<nb_devices;rf++) {
      
      for(unsigned int w = 0; w<nb_wards; w++) {
        
        if(eq_quantity=="variable") {
          if(refill_counter(rf,w)>= refill_freq(rf,w)) {N_mat(rf,w)=refill_quantities(rf,w); refill_counter(rf,w) = 0;}
          else {refill_counter(rf,w)=refill_counter(rf,w)+1;}
          
          if(bin_counter(rf,w)>= equip_bin(rf,w)) {C_mat(rf,w)=0; Npu_mat(rf,w)=0; bin_counter(rf,w) = 0;}
          else {bin_counter(rf,w)=bin_counter(rf,w)+1;}
          
        }
        
        if(refill_counter(rf,w)>= refill_freq(rf,w)){
          
          refill_counter(rf,w) = 0;
          
          if(eq_quantity=="fixed") {
            
            N_mat(rf,w) =  N_mat(rf,w)+refill_quantities(rf,w);
            C_mat(rf,w)=0;
            Npu_mat(rf,w)=0;
            
          }
        }
        
        else {refill_counter(rf,w)=refill_counter(rf,w)+1;}
        
      }
    }
  }
  
  if(output=="simple") {
    
    List L = List::create(_["incidence"] = incidence,
                          _["new_patients"] = new_patient_vec,
                          _["s_patients"] = s_patients,
                          _["cont_new_patients"] = cont_new_patient_vec,
                          _["ward_event"] = ward_event,
                          _["ward_s"] = ward_s,
                          _["count_patients"] = count_patients,
                          _["cont_mat"] = cont_mat,
                          _["used_mat"] = used_mat,
                          _["list_proc_equip"] = list_proc_equip,
                          _["vec_nb_proc"] = vec_nb_proc,
                          _["inf_eq"] = inf_eq,
                          _["count_patient_ward"] = count_patient_ward,
                          _["count_tests"]= count_tests,
                          _["eq_usage"]= eq_usage
                            
    return L;
                          

  }
}

//This function is essentially the same as the one above but allows a "perfect" comparison of intervention vs. baseline scenarios

List BlooPaTH_model_inter(int t, //simulation time
                          int n_patients, // number of patients
                          int nb_wards, // number of wards
                          int nb_adm, // number of admission routes
                          arma::drowvec prev_init, // initial prevalence for each admission route
                          std::string prev_type, //"ward" level or "hospital" level
                          arma::dmat WT_matrix, //transition matrix, specific format (merged matrix for all admission routes) but any time-step will do
                          arma::dmat init_prob, // matrix of probabilities of admission for first ward for both adm routes, nb_adm rows .
                          arma::dcolvec adm_prob,// probabilities of admission for each route
                          int nb_procedures, //number of procedures
                          int nb_equipments, // number of equipments
                          arma::dmat PPM_matrix, // probabilities of undergoing procedures (columns) in each ward (row), merged for both admission routes; Any time-step is ok but need to be the same as the transition matrix
                          Rcpp::DataFrame dist_risk, // distribution parameters of the risk for each procedure (dist can be normal, log-normal or pert)
                          arma::imat nb_material_new, // quantity of equipment for each type in  each ward
                          arma::imat nb_material_used, // quantity of equipment for each type in  each ward
                          arma::imat nb_material_contaminated, // quantity of equipment for each type in  each ward
                          int min_e_phase, // min value for the eclipse phase
                          int max_e_phase, // max value for the eclipse phase
                          float time_step, // must be expressed in hours (ex: if time step is 5min then inform 1/12)
                          bool patient_followup, // useless
                          arma::imat table_proc_equip, // correspondence table between procedures (col 1) and equipments (col 2)
                          StringVector procedure_names, // vector of procedure names
                          StringVector equipment_names, // vector of equipments names
                          arma::drowvec sterilization_prob, // probability of good sterilization for each equipment
                          std::string eq_quantity, // should the quantity be fixed ? "fixed" (then refill quantities not used) or "variable"
                          arma::imat refill_quantities, // quantity of equipment to add to the pool of unused equipment for each type at each refill event
                          arma::fmat refill_freq,// refill frequency for each equipment type (must be expressed in hours)
                          arma::fmat equip_bin,// frequence alaquelle on jete le materiel
                          std::string output="simple",//if output == "all" then returns all tables, if "simple" returns nb of new cases each day and new entry in the hospital
                          float prob_screening=0, // probability of systematic screening upon admission
                          std::string intervention = "none",//none, patient-based or ward-based
                          arma::irowvec at_risk_wards=0, // vector giving the most at risk wards
                          arma::frowvec threshold_reuse=0, // threshold before equipment reuse 
                          std::string pathogen = "HCV",
                          int id_sim=1)
{
  
  List list_proc_equip = transform_list(nb_procedures,table_proc_equip);
  
  
  //clean freq en vecteur pour associer a chaque type de materiel
  arma::dmat storage_adm_TM;
  arma::dmat storage_adm_PPM;
  arma::drowvec storage_init_prob;
  
  StringVector subset_col_proc = table_proc_equip[0];
  StringVector subset_col_equip = table_proc_equip[0];
  
  //Vector of indices to recover each individual transition matrix as the entered matrix is of size ((nb_wards+1)*nb_adm)*(nb_wards+1)
  arma::irowvec index_TM = arma::linspace<arma::irowvec>(0,(nb_wards+1)*nb_adm,nb_adm+1);
  
  //Vector of indices to recover each individual transition matrix as the entered matrix is of size ((nb_wards+1)*nb_adm)*(nb_wards+1)
  arma::irowvec index_PPM = arma::linspace<arma::irowvec>(0,(nb_wards)*nb_adm,nb_adm+1);
  
  //WT_matrix List = different transition matrices associated with different admission places
  // WT_matrix must be same dim than adm_prob
  //adm_prob different from init_prob,adm_prob = admission place and init_prob = first place where proc is performed
  //proc_names, PPM_matrix, and dist risk must be in same order
  //dist_risk = DF with 4 cols, first 3 cols = distribution parameters, fourth col = dist type = "norm", "lnorm" or "pert"
  //PPM_matrix = matrix with nb_proc rows and nb_wards columns : prob of undergoing different kinds of procedures depending on wards
  //gamma = E -> I
  //rho_vec = vector of probs of contaminating material, depending on proc type
  //clean_freq =  frequency at which material is sterilized/thrown away
  //clean_eff = efficiency of material cleaning
  //proc_duration and time step must be in hours
  
  if(prev_type=="ward" && prev_init.size() != nb_wards) {stop("Error: Initial prevalence type chosen is at a ward-level, should be the same size as the number of wards");}
  if(prev_type=="hospital" && prev_init.size() != nb_adm) {stop("Error: Initial prevalence type chosen is at a hospital-level, should be the same size as the number of admission routes");}
  //if(WT_matrix.n_rows < n_wards) {stop("Error: number of rows in WT_matrix < n_wards");}
  
  
  
  //vector for 2= positive, 0=negative
  arma::drowvec pos_neg_vec = {0,2};
  arma::dcolvec cont_neg_vec = {0,1};
  arma::dcolvec failure_success = {0,1};
  arma::drowvec screened_pos_neg_vec = {1,2};
  
  
  //number of patients in the hospital
  int number_pat = n_patients;
  
  //Matrix for status of patients in the hospital
  arma::dmat pop_hosp_status = arma::dmat(n_patients,t+1,arma::fill::zeros); //No need to t+1 because first index = 0
  arma::dmat pop_hosp_status_inter = arma::dmat(n_patients,t+1,arma::fill::zeros); //No need to t+1 because first index = 0
  
  //Matrix for available material through time
  
  arma::icube available_mat_time(nb_equipments,nb_wards,t+1);
  arma::icube available_mat_time_inter(nb_equipments,nb_wards,t+1);
  
  arma::icube cont_mat_time(nb_equipments,nb_wards,t+1);
  arma::icube cont_mat_time_inter(nb_equipments,nb_wards,t+1);
  
  arma::icube used_mat_time(nb_equipments,nb_wards,t+1);
  arma::icube used_mat_time_inter(nb_equipments,nb_wards,t+1); 
  
  
  // counting nb of patients hospitalized on the periode 
  int count_patients = 0;
  
  // counting nb of tests
  int count_tests = 0;
  
  // counting nb of patients hospitalized on the periode in each ward
  arma::irowvec count_patient_ward(nb_wards,arma::fill::zeros);
  
  // counting nb of contaminated patients hospitalized on the periode in each ward
  arma::irowvec cont_new_patient_ward(nb_wards,arma::fill::zeros);
  
  //Matrix for location of patients in the hospital
  arma::dmat pop_hosp_loc = arma::dmat(n_patients,t+1,arma::fill::zeros); //No need to t+1 because first index = 0
  
  //Matrix of procedures undergone by patients during hosp
  arma::dmat pop_proc = arma::dmat(n_patients,t+1,arma::fill::zeros); //No need to t+1 because first index = 0
  
  //Matrix of number of previously used materials (initialized at 0 for each procedure in each ward)
  arma::imat Npu_mat = nb_material_used;
  arma::imat Npu_mat_inter = nb_material_used;
  
  //Matrix of number of available material (initialized at nb_material because no material has been used at start)
  arma::imat N_mat = nb_material_new;
  arma::imat N_mat_inter = nb_material_new;
  
  //Matrix of number of contaminated materials for each ward and each procedure (C compartment)
  arma::imat C_mat = nb_material_contaminated;
  arma::imat C_mat_inter = nb_material_contaminated;
  
  //Matrix storing all departures of patients
  arma::umat patient_departure; //= arma::dmat(2,1,arma::fill::zeros);
  
  //sequence for wards
  arma::drowvec wards_seq = arma::linspace<arma::drowvec>(1,nb_wards,nb_wards);
  
  //sequence for admission places
  arma::drowvec adm_seq = arma::linspace<arma::drowvec>(1,nb_adm,nb_adm);
  
  //sequence for proc
  arma::drowvec proc_seq = arma::linspace<arma::drowvec>(1,nb_procedures,nb_procedures);
  
  //incidence vector
  arma::drowvec incidence(t,arma::fill::zeros);
  arma::drowvec incidence_inter(t,arma::fill::zeros);
  
  //number of s patients vector
  arma::drowvec s_patients(t,arma::fill::zeros);
  arma::drowvec s_patients_inter(t,arma::fill::zeros);
  
  //new patients
  arma::drowvec new_patient_vec(t,arma::fill::zeros);
  new_patient_vec(0) = n_patients;
  
  //contaminated new patients
  arma::drowvec cont_new_patient_vec(t,arma::fill::zeros);
  arma::drowvec cont_new_patient_vec_inter(t,arma::fill::zeros);
  
  //time counter to assign patient in the E state to the I state
  arma::drowvec state_counter(n_patients,arma::fill::zeros);// starts at 0
  arma::drowvec state_counter_inter(n_patients,arma::fill::zeros);// starts at 0
  
  //time counter to assign patient in the E state to the I state
  arma::imat refill_counter(nb_equipments,nb_wards,arma::fill::zeros);// starts at 0
  arma::imat refill_counter_inter(nb_equipments,nb_wards,arma::fill::zeros);// starts at 0
  
  //time counter to assign patient in the E state to the I state
  arma::imat bin_counter(nb_equipments,nb_wards,arma::fill::zeros);// starts at 0
  arma::imat bin_counter_inter(nb_equipments,nb_wards,arma::fill::zeros);// starts at 0
  
  
  //store equipment usage
  arma::imat eq_usage(nb_equipments,nb_wards,arma::fill::zeros);// starts at 0
  arma::imat eq_usage_inter(nb_equipments,nb_wards,arma::fill::zeros);// starts at 0
  
  //time counter to assign patient in the E state to the I state
  arma::drowvec admission_route(n_patients,arma::fill::zeros);// starts at 0
  //arma::drowvec admission_route_inter(n_patients,arma::fill::zeros);// starts at 0
  
  //screened vector upon admission
  arma::drowvec screened(n_patients,arma::fill::zeros);// starts at 0
  
  arma::dcolvec screen_prob_vec = {1-prob_screening,prob_screening};
  
  //Rcout << screen_prob_vec;
  //Matrix to know in patient has been screened for HCV or HBV or HIV , 0 = not screened, 1 = screened and - , 2 = screened and +
  arma::imat ward_screening = arma::imat(n_patients,nb_wards,arma::fill::zeros); //No need to t+1 because first index = 0
  
  //vector to store nb of infection events
  arma::drowvec ward_event(nb_wards,arma::fill::zeros);// starts at 0
  arma::drowvec ward_event_inter(nb_wards,arma::fill::zeros);// starts at 0
  
  //matrix to store nb of suceptible patients in each ward
  arma::dmat ward_s = arma::dmat(nb_wards,t+1,arma::fill::zeros);
  arma::dmat ward_s_inter = arma::dmat(nb_wards,t+1,arma::fill::zeros);
  
  //vector to store nb of infection events
  arma::drowvec cont_mat(nb_equipments,arma::fill::zeros);// starts at 0
  arma::drowvec cont_mat_inter(nb_equipments,arma::fill::zeros);// starts at 0
  
  //vector to store nb of infection events
  arma::drowvec used_mat(nb_equipments,arma::fill::zeros);// starts at 0
  arma::drowvec used_mat_inter(nb_equipments,arma::fill::zeros);// starts at 0
  
  //vector to store nb of proc by proc
  arma::drowvec vec_nb_proc_inter(nb_procedures,arma::fill::zeros);// starts at 0
  arma::drowvec vec_nb_proc(nb_procedures,arma::fill::zeros);// starts at 0
  
  //vector to store nb of infection per equipment type
  arma::drowvec inf_eq(nb_equipments,arma::fill::zeros);// starts at 0
  arma::drowvec inf_eq_inter(nb_equipments,arma::fill::zeros);// starts at 0
  
  //Matrix to know if a patient has already been into a given ward
  arma::imat mat_patient_ward = arma::imat(n_patients,nb_wards,arma::fill::zeros); //No need to t+1 because first index = 0
  
  arma::irowvec init_patient_wards(nb_wards,arma::fill::zeros);
  
  arma::irowvec fill_screened_s = std::vector<int>(nb_wards,1);
  
  arma::irowvec fill_screened_i = std::vector<int>(nb_wards,2);
  
  
  //
  arma::drowvec e_phase = Rcpp::runif( n_patients, min_e_phase, max_e_phase);
  e_phase = round(e_phase);
  
  arma::drowvec patient_id =  arma::linspace<arma::drowvec>(0,n_patients-1,n_patients);
  set_seed(1*id_sim);
  arma::drowvec rdm_patient  = sample(patient_id,n_patients,false);
  
  // For time t=0 initial status are added to pop_hosp
  
  for(int i = 0; i<n_patients;i++) {
    
    
    //set_seed(id_sim); //on reset les seed avant chaque sample pour les differents patients pour s'assurer de la repro de certains processus
    admission_route(i) = sample(adm_seq,1,false,adm_prob)(0);
    //Rcout<<admission_ward(i);
    storage_init_prob= init_prob.row(admission_route(i)-1);
    //Rcout<<wards_seq;
    
    pop_hosp_loc(i,0) = sample(wards_seq,1,false,storage_init_prob.t())(0);
    //Rcout << "ici";
    
    if(prev_type == "hospital") {
      
      arma::dcolvec pos_neg_prob = {1-prev_init((admission_route(i)-1)),prev_init((admission_route(i)-1))};
      
      //test if patient is infected depending on adm route
      
      pop_hosp_status(i,0) = sample(pos_neg_vec,1,false,pos_neg_prob)(0);
      pop_hosp_status_inter(i,0) = sample(pos_neg_vec,1,false,pos_neg_prob)(0);
      
      if(pop_hosp_status(i,0)==2) {cont_new_patient_vec(0) = cont_new_patient_vec(0)+1;}
      if(pop_hosp_status_inter(i,0)==2) {cont_new_patient_vec_inter(0) = cont_new_patient_vec_inter(0)+1;}
      
    }
    
    if(prev_type == "ward") {
      
      arma::dcolvec pos_neg_prob = {1-prev_init((pop_hosp_loc(i,0)-1)),prev_init((pop_hosp_loc(i,0)-1))};
      
      //Rcout<< pos_neg_prob;
      //test if patient is infected depending on adm route
      
      pop_hosp_status(i,0) = sample(pos_neg_vec,1,false,pos_neg_prob)(0);
      pop_hosp_status_inter(i,0) = pop_hosp_status(i,0);
      
      
      if(pop_hosp_status(i,0)==2) {cont_new_patient_vec(0) = cont_new_patient_vec(0)+1;}
      if(pop_hosp_status_inter(i,0)==2) {cont_new_patient_vec_inter(0) = cont_new_patient_vec_inter(0)+1;}
    }
    
    
    
    if(intervention == "ward-based" && std::find(at_risk_wards.begin(), at_risk_wards.end(), (pop_hosp_loc(i,0))) != at_risk_wards.end()) {
      
      int screen = sample(failure_success,1,false,screen_prob_vec)(0);
      
      if (screen==1 && pop_hosp_status_inter(i,0)==0){screened(i)=1;  count_tests = count_tests +1 ;}
      if (screen==1 && pop_hosp_status_inter(i,0)==2){screened(i)=2;  count_tests = count_tests +1 ;}
      
    }
    
    if(intervention =="patient-based") {
      
      //Rcout<< screen_prob_vec;
      
      int screen = sample(failure_success,1,false,screen_prob_vec)(0);
      
      if (screen==1 && pop_hosp_status_inter(i,0)==0) {screened(i)=1; count_tests = count_tests +1 ;}
      if (screen==1 && pop_hosp_status_inter(i,0)==2) {screened(i)=2; count_tests = count_tests +1 ;}
      
    }
    
    
  }
  
  
  
  //adding state for out
  int state_out = nb_wards+1;
  arma::drowvec wards_out_seq = wards_seq;
  wards_out_seq.resize(state_out);
  wards_out_seq(wards_out_seq.size()-1) = state_out;
  
  
  
  
  
  
  
  //starts a time=1 because 0 already informed
  
  for(unsigned int time = 0; time<t;time++) {
    
    // set_seed((time+1)*id_sim); 
    rdm_patient  = sample(patient_id,n_patients,false);
    
    
    for(unsigned int j = 0; j<n_patients;j++) {
      
      
      
      //Rcout<< j;
      
      unsigned int p = rdm_patient(j);
      
      int current_ward = pop_hosp_loc(p,time);
      
      // giving submatrix : rows in span  index_TM(0) --> index_TM(1)-1
      storage_adm_TM = WT_matrix(arma::span(index_TM(admission_route(p)-1),index_TM(admission_route(p))-1),arma::span(0,nb_wards));
      storage_adm_TM =  storage_adm_TM.t();
      arma::dcolvec prob_ward =storage_adm_TM.col(current_ward-1);//minus 1 because first index = 0
      
      // Rcout<< index_PPM(admission_ward(p)-1);
      // Rcout<< index_PPM(admission_ward(p))-1;
      
      storage_adm_PPM = PPM_matrix(arma::span(index_PPM(admission_route(p)-1),index_PPM(admission_route(p))-1),arma::span(0,(nb_procedures-1)));
      storage_adm_PPM =storage_adm_PPM.t();
      //Rcout<< storage_adm_PPM;
      //Rcout << prob_ward;
      
      //set_seed(((p+1)+(time+1))*id_sim); 
      pop_proc(p,time) = sample(proc_seq,1,false, storage_adm_PPM.col(pop_hosp_loc(p,time)-1))(0); // random draw of proc in prob linked to selected ward
      
      //Rcout<<PPM_matrix.col(pop_hosp_loc(p,time)-1);
      
      int current_ward = pop_hosp_loc(p,time);
      
      int current_proc = pop_proc(p,time) ;
      
      int viremic_status = pop_hosp_status(p,time) ; // viremic status for given a patient
      int viremic_status_inter = pop_hosp_status_inter(p,time) ; // viremic status for given a patient for the intervention assessment
      
      
      //count new patient entering ward
      if(mat_patient_ward(p,current_ward-1) == 0) {
        
        count_patient_ward(current_ward-1) = count_patient_ward(current_ward-1)+1;
        
      }
      
      //notifies that patient p has been in ward current_ward
      mat_patient_ward(p,current_ward-1) = 1;
      
      
      if (current_proc != (nb_procedures+1)) {vec_nb_proc(current_proc-1) =vec_nb_proc(current_proc-1)+1;}
      
      //Rcout << ward_screening(p,(current_ward-1));
      
      // if(time>0) {
      //   if(intervention == "ward-based" && 
      //      ward_screening(p,current_ward-1) == 0 && 
      //      std::find(at_risk_wards.begin(), at_risk_wards.end(), current_ward) != at_risk_wards.end() &&
      //      pop_hosp_loc(p,time-1) != pop_hosp_loc(p,time)) {
      //     
      //     int screen = sample(failure_success,1,false,screen_prob_vec)(0);
      //     
      //     if (screen==1 && pop_hosp_status(p,current_ward-1)==0) {ward_screening(p,(current_ward-1))=1;  count_tests = count_tests +1 ;}
      //     if (screen==1 && pop_hosp_status(p,current_ward-1)==2) {ward_screening(p,(current_ward-1))=2;  count_tests = count_tests +1 ;}
      //     
      //   } 
      //   
      // }
      
      
      
      if(time>0) {
        if((intervention == "ward-based") && 
           screened(p) == 0 && 
           std::find(at_risk_wards.begin(), at_risk_wards.end(), current_ward) != at_risk_wards.end() &&
           pop_hosp_loc(p,time-1) != pop_hosp_loc(p,time)) {
          
          int screen = sample(failure_success,1,false,screen_prob_vec)(0);
          
          Rcout<< 'ici';
          
          if (screen==1 && pop_hosp_status_inter(p,time-1)==0) {screened(p) = 1;  count_tests = count_tests +1 ;}
          if (screen==1 && pop_hosp_status_inter(p,time-1)==2) {screened(p) = 2;  count_tests = count_tests +1 ;}
          
        } 
        
      }
      
      //////////////////////////
      // NO INTERVENTION PART //
      //////////////////////////
      
      //Rcout << "ici";
      //for infected patients
      //current_proc != nb_proc because nb_procth procedure = no procedure
      if((viremic_status == 2 or viremic_status == 3) && current_proc != nb_procedures && current_ward != state_out) {
        
        
        //list of involved equipments in given procedure
        arma::irowvec list_equip = list_proc_equip(current_proc-1);
        
        //loop over involved equipments
        for(int eq = 0; eq< list_equip.n_cols; eq++ ) {
          
          
          
          int current_device = list_equip(eq); 
          
          eq_usage(current_device-1,current_ward-1) = eq_usage(current_device-1,current_ward-1) + 1;
          
          float n_available = N_mat(current_device-1,current_ward-1); //nb of total available material for a given proc and ward
          
          float n_pu = Npu_mat(current_device-1,current_ward-1); //nb of total previously used for a given proc and ward
          
          
          //random draw for reused material
          int reuse = 0;
          
          
          // Rcout << "ici 2";
          
          float init_eq = n_available + n_pu;
          
          //if no more available material
          if(n_available < 1) {reuse = 1;}
          
          //if no previously used material
          if(n_pu == 0) {reuse = 0;}
          
          //probability of reused material being contaminated
          
          //if reusing draw to know if reused mat is contaminated
          if(reuse == 1) {
            
            used_mat(current_device-1) = used_mat(current_device-1)+1;
            
            float n_c = C_mat(current_device-1,current_ward-1);
            double prob_Cpw = n_c / n_pu; // prob of drawing a contaminated equipment
            
            // Rcout<< prob_Cpw; 
            if(prob_Cpw > 1) {Rcout<< prob_Cpw; prob_Cpw=1;}
            
            //Rcout<<prob_Cpw ;
            arma::dcolvec vec_prob_Cpw = {1-prob_Cpw,prob_Cpw};
            
            set_seed(((p+1)+(time+1))*id_sim);
            
            int draw_mat = sample(failure_success,1,false,vec_prob_Cpw)(0);
            
            //if material not previously contaminated, it becomes contaminated with a probability equal to sterilization efficiency for the given equipment
            if(draw_mat==0){
              
              // probability of equipment being well sterilized
              double prob_ster = sterilization_prob(current_device-1);
              arma::dcolvec vec_prob_ster = {1-prob_ster,prob_ster};
              int sterilized = sample(failure_success,1,false,vec_prob_ster)(0);
              
              //if patient is detected as positive
              // if sterilization fails
              if(sterilized==0) {C_mat(current_device-1,current_ward-1) = C_mat(current_device-1,current_ward-1)+1;}
              
            }
            
            //if material previously contaminated, it can be sterilized
            if(draw_mat==1){
              
              cont_mat(current_device-1) = cont_mat(current_device-1)+1;
              
              // probability of equipment being well sterilized
              double prob_ster = sterilization_prob(current_device-1);
              arma::dcolvec vec_prob_ster = {1-prob_ster,prob_ster};
              set_seed(((p+1)+(time+1))*id_sim);
              int sterilized = sample(failure_success,1,false,vec_prob_ster)(0);
              
              //if patient is detected as positive
              //if(screened(p) == 2) {sterilized = 1;}
              
              // if sterilization is effective
              if(sterilized==1) {C_mat(current_device-1,current_ward-1) = C_mat(current_device-1,current_ward-1)-1;}
              
            }
            
          }
          
          //if not reusing material then it gets contaminated depending on the good sterilization prob
          if(reuse == 0) {
            
            //adding 1 to matrix of previously used material
            Npu_mat(current_device-1,current_ward-1) = Npu_mat(current_device-1,current_ward-1) + 1;
            
            //removing 1 to matrix of available material
            N_mat(current_device-1,current_ward-1) = N_mat(current_device-1,current_ward-1)-1;
            
            double prob_ster = sterilization_prob(current_device-1);
            arma::dcolvec vec_prob_ster = {1-prob_ster,prob_ster};
            set_seed(((p+1)+(time+1))*id_sim);
            int sterilized = sample(failure_success,1,false,vec_prob_ster)(0);
            
            
            //material becomes contaminated
            if(sterilized==0) {C_mat(current_device-1,current_ward-1) = C_mat(current_device-1,current_ward-1)+1;}
            
          }
          
        }
        
      }
      
      //patients whom are infected but not yet infectious
      if (viremic_status == 1 && current_proc != nb_procedures && current_ward != state_out) {
        
        //list of involved equipments in given procedure
        arma::irowvec list_equip = list_proc_equip(current_proc-1);
        
        //loop over involved equipments
        for(int eq = 0; eq< list_equip.n_cols; eq++ ) {
          
          int current_device = list_equip(eq);
          
          eq_usage(current_device-1,current_ward-1) = eq_usage(current_device-1,current_ward-1) + 1;
          
          float n_available = N_mat(current_device-1,current_ward-1); //nb of total available material for a given proc and ward
          
          float n_pu = Npu_mat(current_device-1,current_ward-1); //nb of total previously used for a given proc and ward
          
          
          //random draw for reused material
          int reuse = 0;
          
          
          float init_eq = n_available + n_pu;
          //Rcout << init_eq;
          
          //if no more available material
          if(n_available < 1) {reuse = 1;}
          
          
          //if no previously used material
          if(n_pu == 0) {reuse = 0;}
          
          //if no reuse
          if(reuse == 0) {
            
            Npu_mat(current_device-1,current_ward-1) = Npu_mat(current_device-1,current_ward-1) + 1;
            N_mat(current_device-1,current_ward-1) = N_mat(current_device-1,current_ward-1)-1;
            
          }
          
          //if reusing draw to know if reused mat is contaminated
          if(reuse == 1) {
            
            
            used_mat(current_device-1) = used_mat(current_device-1)+1;
            
            float n_c = C_mat(current_device-1,current_ward-1);
            double prob_Cpw = n_c / n_pu; // prob of drawing a contaminated equipment
            
            // Rcout<< prob_Cpw;
            
            if(prob_Cpw > 1) {Rcout<< prob_Cpw; prob_Cpw=1;}
            
            //Rcout<<prob_Cpw ;
            arma::dcolvec vec_prob_Cpw = {1-prob_Cpw,prob_Cpw};
            
            set_seed(((p+1)+(time+1))*id_sim);
            
            int draw_mat = sample(failure_success,1,false,vec_prob_Cpw)(0);
            
            //if material contaminated
            if(draw_mat==1){
              
              cont_mat(current_device-1) = cont_mat(current_device-1)+1;
              
              // probability of equipment being well sterilized
              double prob_ster = sterilization_prob(current_device-1);
              arma::dcolvec vec_prob_ster = {1-prob_ster,prob_ster};
              set_seed(((p+1)+(time+1))*id_sim);
              int sterilized = sample(failure_success,1,false,vec_prob_ster)(0);
              
              // if sterilization is well performed
              if(sterilized==1) {C_mat(current_device-1,current_ward-1) = C_mat(current_device-1,current_ward-1)-1;}
              
            }
          }
          
        }
        
        
        
      }
      
      
      //susceptible patients
      if (viremic_status == 0 && current_proc != nb_procedures && current_ward != state_out) {
        
        //list of involved equipments in given procedure
        arma::irowvec list_equip = list_proc_equip(current_proc-1);
        
        // Vector of the prob of equipment being contaminated (used if multiple equipments for one procedure)
        arma::drowvec prob_eq_vec(list_equip.n_cols);
        
        //loop over involved equipments
        for(int eq = 0; eq< list_equip.n_cols; eq++ ) {
          
          int current_device = list_equip(eq); 
          
          eq_usage(current_device-1,current_ward-1) = eq_usage(current_device-1,current_ward-1) + 1;
          
          float n_available = N_mat(current_device-1,current_ward-1); //nb of total available material for a given proc and ward
          
          float n_pu = Npu_mat(current_device-1,current_ward-1); //nb of total previously used for a given proc and ward
          
          //random draw for reused material
          int reuse = 0;
          
          
          //if no more available material
          float init_eq = n_available + n_pu;
          //Rcout << init_eq;
          
          if(n_available < 1) {reuse = 1;}
          
          if(n_available > 1 && n_available/init_eq<threshold_reuse(current_device-1)) {
            
            double prob_reuse= 1-(n_available/init_eq);
            
            arma::dcolvec vec_prob_reuse = {1-prob_reuse,prob_reuse};
            
            //Rcout << vec_prob_reuse;
            
            set_seed(((p+1)+(time+1))*id_sim);
            
            reuse = sample(failure_success,1,false,vec_prob_reuse)(0);}
          
          //if no previously used material
          if(n_pu == 0) {reuse = 0;}
          
          
          //probability of reused material being contaminated
          
          //if reusing draw to know if reused mat is contaminated
          if(reuse == 1) {
            
            
            
            used_mat(current_device-1) = used_mat(current_device-1)+1;
            
            
            float n_c = C_mat(current_device-1,current_ward-1);
            double prob_Cpw = n_c / n_pu; // prob of drawing a contaminated equipment
            
            //Rcout<< prob_Cpw;
            
            if(prob_Cpw > 1) {Rcout<< prob_Cpw; prob_Cpw=1;}
            
            
            
            //Rcout<<prob_Cpw ;
            
            arma::dcolvec vec_prob_Cpw = {1-prob_Cpw,prob_Cpw};
            
            set_seed(((p+1)+(time+1))*id_sim);
            
            int draw_mat = sample(failure_success,1,false,vec_prob_Cpw)(0);
            
            prob_eq_vec(eq) = draw_mat;
            
            //if material contaminated
            if(draw_mat==1){
              
              cont_mat(current_device-1) = cont_mat(current_device-1)+1;
              
              // probability of equipment being well sterilized
              double prob_ster = sterilization_prob(current_device-1);
              arma::dcolvec vec_prob_ster = {1-prob_ster,prob_ster};
              set_seed(((p+1)+(time+1))*id_sim);
              int sterilized = sample(failure_success,1,false,vec_prob_ster)(0);
              
              // if sterilization is as success
              if(sterilized==1) {C_mat(current_device-1,current_ward-1) = C_mat(current_device-1,current_ward-1)-1;}
              
            }
            
          }
          
          //no reuse
          if(reuse == 0) {
            
            Npu_mat(current_device-1,current_ward-1) = Npu_mat(current_device-1,current_ward-1) + 1;
            N_mat(current_device-1,current_ward-1) = N_mat(current_device-1,current_ward-1)-1;
          }
          
        }
        
        
        //Rcpp:Rcout << time;
        
        //Rcout<< prob_eq_vec;
        
        //if a least 1 eq is contaminated
        if(sum(prob_eq_vec) >= 1) {
          
          
          StringVector vec_dist = dist_risk["dist"];
          
          if(vec_dist(current_proc-1)=="lnorm") {
            
            
            arma::fvec lmean_vec = dist_risk[0];
            
            arma::fvec lsd_vec = dist_risk[1];
            
            float lmean_risk = lmean_vec(current_proc-1); //parameter 1 for lognormal dist
            
            float lsd_risk = lsd_vec(current_proc-1); //parameter 2 for lognormal dist
            
            float risk_p = Rcpp::as<float> (random_lnorm(_["n"]=1,_["meanlog"]=lmean_risk,_["sdlog"]=lsd_risk)); //risk of getting HCV contaminated if material is contaminated
            
            float prob_E =risk_p ; 
            
            if(pathogen =="HBV") {
              
              float HBV_mult = Rcpp::as<float> (random_unif(_["n"]=1,_["min"]=4,_["max"]=20)); //risk mult for HBV
              
              prob_E = risk_p*HBV_mult;
              
              Rcout<< prob_E;
              
            }
            
            
            
            Rcout<< prob_E;
            
            //Probability of going from the S to the E state
            
            
            //float prob_E = risk_p;
            
            if (prob_E<0) {prob_E=0;}
            if (prob_E>1) {prob_E=1;}
            
            
            arma::dcolvec prob_inf_vector = {1-prob_E,prob_E};
            
            set_seed(((p+1)+(time+1))*id_sim);
            
            int draw_infection = sample(cont_neg_vec,1,false,prob_inf_vector)(0);
            
            if (draw_infection == 1) {pop_hosp_status(p,time) = 1; incidence(time)=incidence(time)+1; ward_event(current_ward-1)= ward_event(current_ward-1)+1; 
            
            for(int ei = 0; ei < list_equip.n_cols; ei++) {if(prob_eq_vec(ei)==1) {inf_eq(list_equip(ei)-1) = inf_eq(list_equip(ei)-1)+1;} } 
            
            }
            
          }
          
          if(vec_dist(current_proc-1)=="norm") {
            
            arma::fvec mean_vec = dist_risk[0];
            
            arma::fvec sd_vec = dist_risk[1];
            
            float mean_risk = mean_vec(current_proc-1); //parameter 1 for normal dist
            
            float sd_risk = sd_vec(current_proc-1); //parameter 2 for normal dist
            
            float risk_p = Rcpp::as<float> (random_norm(_["n"]=1,_["mean"]=mean_risk,_["sd"]=sd_risk)); //risk of getting HCV contaminated if material is contaminated
            
            float prob_E = risk_p;
            
            if(pathogen=="HBV") {
              
              
              
              float HBV_mult = Rcpp::as<float> (random_unif(_["n"]=1,_["min"]=4,_["max"]=20)); //risk mult for HBV
              
              prob_E = risk_p*HBV_mult;
              
              
              
            }
            
            //float prob_E = risk_p;
            
            if (prob_E<0) {prob_E=0;}
            if (prob_E>1) {prob_E=1;}
            
            
            
            arma::dcolvec prob_inf_vector = {1-prob_E,prob_E};
            
            set_seed(((p+1)+(time+1))*id_sim);
            
            int draw_infection = sample(cont_neg_vec,1,false,prob_inf_vector)(0);
            
            if (draw_infection == 1) {pop_hosp_status(p,time) = 1; incidence(time)=incidence(time)+1; ward_event(current_ward-1)= ward_event(current_ward-1)+1;
            
            for(int ei = 0; ei < list_equip.n_cols; ei++) {if(prob_eq_vec(ei)==1) {inf_eq(list_equip(ei)) = inf_eq(list_equip(ei)-1)+1;} } 
            }
            
          }
          
          if(vec_dist(current_proc-1)=="pert") {
            
            arma::fvec min_vec = dist_risk[0];
            
            arma::fvec mode_vec = dist_risk[1];
            
            arma::fvec max_vec = dist_risk[2];
            
            float min_risk = min_vec(current_proc-1); //parameter 1 for pert dist
            
            float mode_risk = mode_vec(current_proc-1); //parameter 2 for pert dist
            
            float max_risk = max_vec(current_proc-1); //parameter 3 for pert dist
            
            float risk_p = Rcpp::as<float> (random_pert(_["n"]=1,_["min"]=min_risk,_["mode"]=mode_risk,_["max"]=max_risk))/100; //risk of getting HCV contaminated if material is contaminated
            
            float prob_E = risk_p;
            
            if(pathogen =="HBV") {
              
              float HBV_mult = Rcpp::as<float> (random_unif(_["n"]=1,_["min"]=4,_["max"]=20)); //risk mult for HBV
              
              prob_E = risk_p*HBV_mult;
              
            }
            //float prob_E = risk_p;
            
            if (prob_E<0) {prob_E=0;}
            if (prob_E>1) {prob_E=1;}
            
            
            
            arma::dcolvec prob_inf_vector = {1-prob_E,prob_E};
            
            set_seed(((p+1)+(time+1))*id_sim);
            
            int draw_infection = sample(cont_neg_vec,1,false,prob_inf_vector)(0);
            
            if (draw_infection == 1) {pop_hosp_status(p,time) = 1; incidence(time)=incidence(time)+1; ward_event(current_ward-1)= ward_event(current_ward-1)+1;
            
            for(int ei = 0; ei < list_equip.n_cols; ei++) {if(prob_eq_vec(ei)==1) {inf_eq(list_equip(ei)-1) = inf_eq(list_equip(ei)-1)+1;} } 
            }
            
          }
          
          
        }
        
      }
      
      
      ////////////////////////////
      // WITH INTERVENTION PART //
      ////////////////////////////
      
      //Rcout << "ici";
      //for infected patients
      //current_proc != nb_proc because nb_procth procedure = no procedure
      if((viremic_status_inter == 2 or viremic_status_inter == 3) && current_proc != nb_procedures && current_ward != state_out) {
        
        
        //list of involved equipments in given procedure
        arma::irowvec list_equip = list_proc_equip(current_proc-1);
        
        //loop over involved equipments
        for(int eq = 0; eq< list_equip.n_cols; eq++ ) {
          
          
          int current_device = list_equip(eq); 
          
          eq_usage_inter(current_device-1,current_ward-1) = eq_usage_inter(current_device-1,current_ward-1) + 1;
          
          float n_available = N_mat_inter(current_device-1,current_ward-1); //nb of total available material for a given proc and ward
          
          float n_pu = Npu_mat_inter(current_device-1,current_ward-1); //nb of total previously used for a given proc and ward
          
          
          //random draw for reused material
          int reuse = 0;
          
          
          // Rcout << "ici 2";
          
          float init_eq = n_available + n_pu;
          
          //if no more available material
          if(n_available < 1) {reuse = 1;}
          
          //if no previously used material
          if(n_pu == 0) {reuse = 0;}
          
          //probabiliused_matty of reused material being contaminated
          
          //if reusing draw to know if reused mat is contaminated
          if(reuse == 1) {
            
            used_mat_inter(current_device-1) = used_mat_inter(current_device-1)+1;
            
            float n_c = C_mat_inter(current_device-1,current_ward-1);
            double prob_Cpw = n_c / n_pu; // prob of drawing a contaminated equipment
            
            // Rcout<< prob_Cpw; 
            
            //Rcout<<prob_Cpw ;
            arma::dcolvec vec_prob_Cpw = {1-prob_Cpw,prob_Cpw};
            
            set_seed(((p+1)+(time+1))*id_sim);
            
            int draw_mat = sample(failure_success,1,false,vec_prob_Cpw)(0);
            
            //if material not previously contaminated, it becomes contaminated with a probability equal to sterilization efficiency for the given equipment
            if(draw_mat==0){
              
              // probability of equipment being well sterilized
              double prob_ster = sterilization_prob(current_device-1);
              arma::dcolvec vec_prob_ster = {1-prob_ster,prob_ster};
              set_seed(((p+1)+(time+1))*id_sim);
              int sterilized = sample(failure_success,1,false,vec_prob_ster)(0);
              
              //if patient is detected as positive
              if( (intervention == "patient-based" or intervention == "ward-based") && screened(p) == 2) {sterilized = 1;}
              
              // if sterilization fails
              if(sterilized==0) {C_mat_inter(current_device-1,current_ward-1) = C_mat_inter(current_device-1,current_ward-1)+1;}
              
            }
            
            //if material previously contaminated, it can be sterilized
            if(draw_mat==1){
              
              cont_mat_inter(current_device-1) = cont_mat_inter(current_device-1)+1;
              
              // probability of equipment being well sterilized
              double prob_ster = sterilization_prob(current_device-1);
              arma::dcolvec vec_prob_ster = {1-prob_ster,prob_ster};
              set_seed(((p+1)+(time+1))*id_sim);
              int sterilized = sample(failure_success,1,false,vec_prob_ster)(0);
              
              if( (intervention == "patient-based" or intervention == "ward-based") && screened(p) == 2) {sterilized = 1;}
              
              //if patient is detected as positive
              //if(screened(p) == 2) {sterilized = 1;}
              
              // if sterilization is effective
              if(sterilized==1) {C_mat_inter(current_device-1,current_ward-1) = C_mat_inter(current_device-1,current_ward-1)-1;}
              
            }
            
          }
          
          //if not reusing material then it gets contaminated depending on the good sterilization prob
          if(reuse == 0) {
            
            //adding 1 to matrix of previously used material
            Npu_mat_inter(current_device-1,current_ward-1) = Npu_mat_inter(current_device-1,current_ward-1) + 1;
            
            //removing 1 to matrix of available material
            N_mat_inter(current_device-1,current_ward-1) = N_mat_inter(current_device-1,current_ward-1)-1;
            
            double prob_ster = sterilization_prob(current_device-1);
            arma::dcolvec vec_prob_ster = {1-prob_ster,prob_ster};
            set_seed(((p+1)+(time+1))*id_sim);
            int sterilized = sample(failure_success,1,false,vec_prob_ster)(0);
            
            if( (intervention == "patient-based" or intervention == "ward-based") && screened(p) == 2) {sterilized = 1;}
            
            
            //material becomes contaminated
            if(sterilized==0) {C_mat_inter(current_device-1,current_ward-1) = C_mat_inter(current_device-1,current_ward-1)+1;}
            
          }
          
        }
        
      }
      
      //patients that are infected but not infectious
      if (viremic_status_inter == 1 && current_proc != nb_procedures && current_ward != state_out) {
        
        //list of involved equipments in given procedure
        arma::irowvec list_equip = list_proc_equip(current_proc-1);
        
        //loop over involved equipments
        for(int eq = 0; eq< list_equip.n_cols; eq++ ) {
          
          int current_device = list_equip(eq);
          
          eq_usage_inter(current_device-1,current_ward-1) = eq_usage_inter(current_device-1,current_ward-1) + 1;
          
          float n_available = N_mat_inter(current_device-1,current_ward-1); //nb of total available material for a given proc and ward
          
          float n_pu = Npu_mat_inter(current_device-1,current_ward-1); //nb of total previously used for a given proc and ward
          
          
          //random draw for reused material
          int reuse = 0;
          
          // if(time>0) {
          // if(intervention == "ward-based" && 
          //    ward_screening(p,current_ward-1) == 0 && 
          //    std::find(at_risk_wards.begin(), at_risk_wards.end(), current_ward) != at_risk_wards.end() &&
          //    pop_hosp_loc(p,time-1) != pop_hosp_loc(p,time)) {
          //   
          //   int screen = sample(failure_success,1,false,screen_prob_vec)(0);
          //   
          //   if (screen==1 && pop_hosp_status(p,current_ward-1)==0) {ward_screening(p,(pop_hosp_loc(p,current_ward-1)-1))=1;  count_tests = count_tests +1 ;}
          //   if (screen==1 && pop_hosp_status(p,current_ward-1)==2) {ward_screening(p,(pop_hosp_loc(p,current_ward-1)-1))=2;  count_tests = count_tests +1 ;}
          //   
          // }
          // 
          // }
          float init_eq = n_available + n_pu;
          //Rcout << init_eq;
          
          //if no more available material
          if(n_available < 1) {reuse = 1;}
          
          //if no previously used material
          if(n_pu == 0) {reuse = 0;}
          
          //if no reuse
          if(reuse == 0) {
            
            Npu_mat_inter(current_device-1,current_ward-1) = Npu_mat_inter(current_device-1,current_ward-1) + 1;
            N_mat_inter(current_device-1,current_ward-1) = N_mat_inter(current_device-1,current_ward-1)-1;
            
          }
          
          //if reusing draw to know if reused mat is contaminated
          if(reuse == 1) {
            
            
            used_mat_inter(current_device-1) = used_mat_inter(current_device-1)+1;
            
            float n_c = C_mat_inter(current_device-1,current_ward-1);
            double prob_Cpw = n_c / n_pu; // prob of drawing a contaminated equipment
            
            // Rcout<< prob_Cpw;
            
            if(prob_Cpw > 1) {Rcout<< prob_Cpw; prob_Cpw=1;}
            
            //Rcout<<prob_Cpw ;
            arma::dcolvec vec_prob_Cpw = {1-prob_Cpw,prob_Cpw};
            
            set_seed(((p+1)+(time+1))*id_sim);
            
            int draw_mat = sample(failure_success,1,false,vec_prob_Cpw)(0);
            
            //if material contaminated
            if(draw_mat==1){
              
              cont_mat_inter(current_device-1) = cont_mat_inter(current_device-1)+1;
              
              // probability of equipment being well sterilized
              double prob_ster = sterilization_prob(current_device-1);
              arma::dcolvec vec_prob_ster = {1-prob_ster,prob_ster};
              set_seed(((p+1)+(time+1))*id_sim);
              int sterilized = sample(failure_success,1,false,vec_prob_ster)(0);
              
              // if sterilization is well performed
              if(sterilized==1) {C_mat_inter(current_device-1,current_ward-1) = C_mat_inter(current_device-1,current_ward-1)-1;}
              
            }
          }
          
        }
        
        
        
      }
      
      
      //susceptible patients
      if (viremic_status_inter == 0 && current_proc != nb_procedures && current_ward != state_out) {
        
        //list of involved equipments in given procedure
        arma::irowvec list_equip = list_proc_equip(current_proc-1);
        
        // Vector of the prob of equipment being contaminated (used if multiple equipments for one procedure)
        arma::drowvec prob_eq_vec(list_equip.n_cols);
        
        //loop over involved equipments
        for(int eq = 0; eq< list_equip.n_cols; eq++ ) {
          
          int current_device = list_equip(eq); 
          
          eq_usage_inter(current_device-1,current_ward-1) = eq_usage_inter(current_device-1,current_ward-1) + 1;
          
          float n_available = N_mat_inter(current_device-1,current_ward-1); //nb of total available material for a given proc and ward
          
          float n_pu = Npu_mat_inter(current_device-1,current_ward-1); //nb of total previously used for a given proc and ward
          
          //random draw for reused material
          int reuse = 0;
          
          // if(time>0) {
          // if(intervention == "ward-based" && 
          //    ward_screening(p,current_ward-1) == 0 && 
          //    std::find(at_risk_wards.begin(), at_risk_wards.end(), current_ward) != at_risk_wards.end() &&
          //    pop_hosp_loc(p,time-1) != pop_hosp_loc(p,time)) {
          //   
          //   int screen = sample(failure_success,1,false,screen_prob_vec)(0);
          //   
          //   if (screen==1 && pop_hosp_status(p,current_ward-1)==0) {ward_screening(p,(pop_hosp_loc(p,current_ward-1)-1))=1;  count_tests = count_tests +1 ;}
          //   if (screen==1 && pop_hosp_status(p,current_ward-1)==2) {ward_screening(p,(pop_hosp_loc(p,current_ward-1)-1))=2;  count_tests = count_tests +1 ;}
          //   
          // }
          // 
          // }
          // 
          //if no more available material
          float init_eq = n_available + n_pu;
          //Rcout << init_eq;
          
          if(n_available < 1) {reuse = 1;}
          
          if(n_available > 1 && n_available/init_eq<threshold_reuse(current_device-1)) {
            
            double prob_reuse= 1-(n_available/init_eq);
            
            arma::dcolvec vec_prob_reuse = {1-prob_reuse,prob_reuse};
            
            //Rcout << vec_prob_reuse;
            
            set_seed(((p+1)+(time+1))*id_sim);
            
            reuse = sample(failure_success,1,false,vec_prob_reuse)(0);}
          
          //if no previously used material
          if(n_pu == 0) {reuse = 0;}
          
          
          //probability of reused material being contaminated
          
          //if reusing draw to know if reused mat is contaminated
          if(reuse == 1) {
            
            
            
            used_mat_inter(current_device-1) = used_mat_inter(current_device-1)+1;
            
            
            float n_c = C_mat_inter(current_device-1,current_ward-1);
            double prob_Cpw = n_c / n_pu; // prob of drawing a contaminated equipment
            
            //Rcout<< prob_Cpw;
            
            if(prob_Cpw > 1) {Rcout<< prob_Cpw; prob_Cpw=1;}
            
            
            
            //Rcout<<prob_Cpw ;
            
            arma::dcolvec vec_prob_Cpw = {1-prob_Cpw,prob_Cpw};
            
            set_seed(((p+1)+(time+1))*id_sim);
            
            int draw_mat = sample(failure_success,1,false,vec_prob_Cpw)(0);
            
            prob_eq_vec(eq) = draw_mat;
            
            //if material contaminated
            if(draw_mat==1){
              
              cont_mat_inter(current_device-1) = cont_mat_inter(current_device-1)+1;
              
              // probability of equipment being well sterilized
              double prob_ster = sterilization_prob(current_device-1);
              arma::dcolvec vec_prob_ster = {1-prob_ster,prob_ster};
              set_seed(((p+1)+(time+1))*id_sim);
              int sterilized = sample(failure_success,1,false,vec_prob_ster)(0);
              
              // if sterilization is as success
              if(sterilized==1) {C_mat_inter(current_device-1,current_ward-1) = C_mat_inter(current_device-1,current_ward-1)-1;}
              
            }
            
          }
          
          //no reuse
          if(reuse == 0) {
            
            Npu_mat_inter(current_device-1,current_ward-1) = Npu_mat_inter(current_device-1,current_ward-1) + 1;
            N_mat_inter(current_device-1,current_ward-1) = N_mat_inter(current_device-1,current_ward-1)-1;
          }
          
        }
        
        
        //Rcpp:Rcout << time;
        
        //Rcout<< prob_eq_vec;
        
        //if a least 1 eq is contaminated
        if(sum(prob_eq_vec) >= 1) {
          
          
          StringVector vec_dist = dist_risk["dist"];
          
          if(vec_dist(current_proc-1)=="lnorm") {
            
            
            arma::fvec lmean_vec = dist_risk[0];
            
            arma::fvec lsd_vec = dist_risk[1];
            
            float lmean_risk = lmean_vec(current_proc-1); //parameter 1 for lognormal dist
            
            float lsd_risk = lsd_vec(current_proc-1); //parameter 2 for lognormal dist
            
            float risk_p = Rcpp::as<float> (random_lnorm(_["n"]=1,_["meanlog"]=lmean_risk,_["sdlog"]=lsd_risk)); //risk of getting HCV contaminated if material is contaminated
            
            float prob_E =risk_p ; 
            
            if(pathogen =="HBV") {
              
              float HBV_mult = Rcpp::as<float> (random_unif(_["n"]=1,_["min"]=4,_["max"]=20)); //risk mult for HBV
              
              prob_E = risk_p*HBV_mult;
              
              Rcout<< prob_E;
              
            }
            
            
            
            Rcout<< prob_E;
            
            //Probability of going from the S to the E state
            
            
            //float prob_E = risk_p;
            
            if (prob_E<0) {prob_E=0;}
            if (prob_E>1) {prob_E=1;}
            
            
            arma::dcolvec prob_inf_vector = {1-prob_E,prob_E};
            
            set_seed(((p+1)+(time+1))*id_sim);
            
            int draw_infection = sample(cont_neg_vec,1,false,prob_inf_vector)(0);
            
            if (draw_infection == 1) {pop_hosp_status_inter(p,time) = 1; incidence_inter(time)=incidence_inter(time)+1; ward_event_inter(current_ward-1)= ward_event_inter(current_ward-1)+1; 
            
            for(int ei = 0; ei < list_equip.n_cols; ei++) {if(prob_eq_vec(ei)==1) {inf_eq_inter(list_equip(ei)-1) = inf_eq_inter(list_equip(ei)-1)+1;} } 
            
            }
            
          }
          
          if(vec_dist(current_proc-1)=="norm") {
            
            arma::fvec mean_vec = dist_risk[0];
            
            arma::fvec sd_vec = dist_risk[1];
            
            float mean_risk = mean_vec(current_proc-1); //parameter 1 for normal dist
            
            float sd_risk = sd_vec(current_proc-1); //parameter 2 for normal dist
            
            float risk_p = Rcpp::as<float> (random_norm(_["n"]=1,_["mean"]=mean_risk,_["sd"]=sd_risk)); //risk of getting HCV contaminated if material is contaminated
            
            float prob_E = risk_p;
            
            if(pathogen=="HBV") {
              
              
              
              float HBV_mult = Rcpp::as<float> (random_unif(_["n"]=1,_["min"]=4,_["max"]=20)); //risk mult for HBV
              
              prob_E = risk_p*HBV_mult;
              
              
              
            }
            
            //float prob_E = risk_p;
            
            if (prob_E<0) {prob_E=0;}
            if (prob_E>1) {prob_E=1;}
            
            
            
            arma::dcolvec prob_inf_vector = {1-prob_E,prob_E};
            
            set_seed(((p+1)+(time+1))*id_sim);
            
            int draw_infection = sample(cont_neg_vec,1,false,prob_inf_vector)(0);
            
            if (draw_infection == 1) {pop_hosp_status_inter(p,time) = 1; incidence_inter(time)=incidence_inter(time)+1; ward_event_inter(current_ward-1)= ward_event_inter(current_ward-1)+1;
            
            for(int ei = 0; ei < list_equip.n_cols; ei++) {if(prob_eq_vec(ei)==1) {inf_eq_inter(list_equip(ei)) = inf_eq_inter(list_equip(ei)-1)+1;} } 
            }
            
          }
          
          if(vec_dist(current_proc-1)=="pert") {
            
            arma::fvec min_vec = dist_risk[0];
            
            arma::fvec mode_vec = dist_risk[1];
            
            arma::fvec max_vec = dist_risk[2];
            
            float min_risk = min_vec(current_proc-1); //parameter 1 for pert dist
            
            float mode_risk = mode_vec(current_proc-1); //parameter 2 for pert dist
            
            float max_risk = max_vec(current_proc-1); //parameter 3 for pert dist
            
            float risk_p = Rcpp::as<float> (random_pert(_["n"]=1,_["min"]=min_risk,_["mode"]=mode_risk,_["max"]=max_risk))/100; //risk of getting HCV contaminated if material is contaminated
            
            float prob_E = risk_p;
            
            if(pathogen =="HBV") {
              
              float HBV_mult = Rcpp::as<float> (random_unif(_["n"]=1,_["min"]=4,_["max"]=20)); //risk mult for HBV
              
              prob_E = risk_p*HBV_mult;
              
            }
            //float prob_E = risk_p;
            
            if (prob_E<0) {prob_E=0;}
            if (prob_E>1) {prob_E=1;}
            
            
            
            arma::dcolvec prob_inf_vector = {1-prob_E,prob_E};
            
            set_seed(((p+1)+(time+1))*id_sim);
            
            int draw_infection = sample(cont_neg_vec,1,false,prob_inf_vector)(0);
            
            if (draw_infection == 1) {pop_hosp_status_inter(p,time) = 1; incidence_inter(time)=incidence_inter(time)+1; ward_event_inter(current_ward-1)= ward_event_inter(current_ward-1)+1;
            
            for(int ei = 0; ei < list_equip.n_cols; ei++) {if(prob_eq_vec(ei)==1) {inf_eq_inter(list_equip(ei)-1) = inf_eq_inter(list_equip(ei)-1)+1;} } 
            }
            
          }
          
          
        }
        
      }
      
      
      
      
      // available_mat_time.col(time+1) = sum(N_mat,1);
      
      // used_mat_time.col(time+1) = sum(Npu_mat,1);
      // contaminated_mat_time.col(time+1) = sum(C_mat,1);
      
      
      //adding 1 to status counter
      if(pop_hosp_status(p,time) == 1) {state_counter(p) = state_counter(p) + 1;}
      if(pop_hosp_status_inter(p,time) == 1) {state_counter_inter(p) = state_counter_inter(p) + 1;}
      
      //if counter exceed the time of the eclipse phase then patient goes from E to I
      if(state_counter(p) >= e_phase(p) && pop_hosp_status(p,time) == 1 ) {pop_hosp_status(p,time) = 3; }
      if(state_counter_inter(p) >= e_phase(p) && pop_hosp_status_inter(p,time) == 1 ) {pop_hosp_status_inter(p,time) = 3; }
      
      //set_seed(((p+1)+(time+1))*id_sim);
      pop_hosp_loc(p,time+1) = sample(wards_out_seq,1,false,prob_ward)(0); //first ward for new patient
      
      
      pop_hosp_status(p,time+1) = pop_hosp_status(p,time); //patient has same status than time before
      pop_hosp_status_inter(p,time+1) = pop_hosp_status_inter(p,time); //patient has same status than time before
      
      //pop_proc(p,time+1) = pop_proc(p,time);
      
      
      if(pop_hosp_status(p,time)==0) {s_patients(time)=s_patients(time)+1; ward_s(current_ward-1,time)=ward_s(current_ward-1,time)+1;}
      if(pop_hosp_status_inter(p,time)==0) {s_patients_inter(time)=s_patients_inter(time)+1; ward_s_inter(current_ward-1,time)=ward_s_inter(current_ward-1,time)+1;}
      
      
      // test_eq(time) = N_mat(4,1);
      // 
      // test_eq_u(time) = Npu_mat(4,1) - C_mat(4,1); 
      // 
      // test_eq_cont(time) =C_mat(4,1);
      
      // if a patient leaves the hospital, another one is entering
      if( pop_hosp_loc(p,time+1) == state_out) {
        
        //initialize to zero so that history of passing through wards is clear for new patient
        mat_patient_ward.row(p) = init_patient_wards;
        
        ward_screening.row(p) = init_patient_wards;
        
        count_patients = count_patients +1;
        
        
        new_patient_vec(time) = new_patient_vec(time)+1;
        
        arma::ucolvec new_departure = {time+1,p};
        patient_departure = join_horiz(patient_departure ,new_departure) ;
        
        //set_seed(((p+1)+(time+1))*id_sim);
        admission_route(p) = sample(adm_seq,1,false,adm_prob)(0);
        storage_init_prob= init_prob.row(admission_route(p)-1);
        pop_hosp_loc(p,time+1) = sample(wards_seq,1,false,storage_init_prob.t())(0);//if patient leaves the hospital, another patient instantly replaces him
        
        //count_patient_ward(pop_hosp_loc(p,time+1)-1) =   count_patient_ward(pop_hosp_loc(p,time+1)-1)+1;      
        
        state_counter(p) = 0; //if patient leaves the hospital, counter of infection for patient p goes back to 0
        state_counter_inter(p) = 0;
        
        if(prev_type == "ward") {
          
          arma::dcolvec pos_neg_prob = {1-prev_init((pop_hosp_loc(p,time+1)-1)),prev_init((pop_hosp_loc(p,time+1)-1))};
          
          //test if patient is infected depending on adm route
          //set_seed(((p+1)+(time+1))*id_sim);
          pop_hosp_status(p,time+1) = sample(pos_neg_vec,1,false,pos_neg_prob)(0);
          pop_hosp_status_inter(p,time+1) = pop_hosp_status(p,time+1);
          
        }
        
        if(prev_type == "hospital") {
          
          arma::dcolvec pos_neg_prob = {1-prev_init((admission_route(p)-1)),prev_init((admission_route(p)-1))};
          
          //test if patient is infected depending on adm route
          //set_seed(((p+1)+(time+1))*id_sim);
          pop_hosp_status(p,time+1) = sample(pos_neg_vec,1,false,pos_neg_prob)(0);
          pop_hosp_status_inter(p,time+1) = pop_hosp_status(p,time+1);
        }
        
        e_phase(p) = round(Rcpp::runif(1,min_e_phase,max_e_phase)(0));
        
        
        screened(p)=0;
        
        if(intervention == "ward-based" && std::find(at_risk_wards.begin(), at_risk_wards.end(), (pop_hosp_loc(p,time+1))) != at_risk_wards.end()) {
          
          int screen = sample(failure_success,1,false,screen_prob_vec)(0);
          
          if (screen==1 && pop_hosp_status_inter(p,time+1)==0) {screened(p)=1; count_tests = count_tests +1 ;}
          if (screen==1 && pop_hosp_status_inter(p,time+1)==2) {screened(p)=2; count_tests = count_tests +1 ;}
          
        }
        
        if(intervention =="patient-based") {
          
          int screen_new = sample(failure_success,1,false,screen_prob_vec)(0);
          
          // Rcout << "ici 3 ";
          
          if (screen_new==1 && pop_hosp_status_inter(p,time+1)==0) {screened(p)=1;  count_tests = count_tests +1 ;}
          
          if (screen_new==1 && pop_hosp_status_inter(p,time+1)==2) {screened(p)=2;  count_tests = count_tests +1 ;}
          
        }
        
      }
      
      
      //count_patient_ward(pop_hosp_loc(p,time+1)-1) =   count_patient_ward(pop_hosp_loc(p,time+1)-1)+1;
      
      
      // if ( std::find(vec.begin(), vec.end(), item) != vec.end() )
      
      // if(pop_hosp_status(p,time)==0 &&  mat_patient_ward(p,current_ward-1) == 0) {
      // 
      //   count_patient_ward(current_ward-1) = count_patient_ward(current_ward-1)+1;
      // 
      // }
      
      // if((pop_hosp_status(p,time)==1 or pop_hosp_status(p,time)==2 or pop_hosp_status(p,time)==3)
      //      &&  mat_patient_ward(p,current_ward-1) == 0) {
      // 
      //   count_patient_ward(current_ward-1) = count_patient_ward(current_ward-1)+1;
      // 
      // }
      
      
    }
    
    
    
    
    
    
    // check if there s a need to refill equipment or throw away equipment
    for(unsigned int rf=0;rf<nb_equipments;rf++) {
      
      for(unsigned int w = 0; w<nb_wards; w++) {
        
        if(eq_quantity=="variable") {
          if(refill_counter(rf,w)>= refill_freq(rf,w)) {N_mat(rf,w)=refill_quantities(rf,w); refill_counter(rf,w) = 0;}
          else {refill_counter(rf,w)=refill_counter(rf,w)+1;}
          
          if(refill_counter_inter(rf,w)>= refill_freq(rf,w)) {N_mat_inter(rf,w)=refill_quantities(rf,w); refill_counter_inter(rf,w) = 0;}
          else {refill_counter_inter(rf,w)=refill_counter_inter(rf,w)+1;}
          
          if(bin_counter(rf,w)>= equip_bin(rf,w)) {C_mat(rf,w)=0; Npu_mat(rf,w)=0; bin_counter(rf,w) = 0;}
          else {bin_counter(rf,w)=bin_counter(rf,w)+1;}
          
          if(bin_counter_inter(rf,w)>= equip_bin(rf,w)) {C_mat_inter(rf,w)=0; Npu_mat_inter(rf,w)=0; bin_counter_inter(rf,w) = 0;}
          else {bin_counter_inter(rf,w)=bin_counter_inter(rf,w)+1;}
          
        }
        
        if(refill_counter(rf,w)>= refill_freq(rf,w)){
          
          refill_counter(rf,w) = 0;
          
          if(eq_quantity=="fixed") {
            
            //N_mat(rf,w) =  N_mat(rf,w) + Npu_mat(rf,w);
            N_mat(rf,w) =  N_mat(rf,w)+refill_quantities(rf,w);
            C_mat(rf,w)=0;
            Npu_mat(rf,w)=0;
            
          }
          
        }
        
        else {refill_counter(rf,w)=refill_counter(rf,w)+1;}
        
        
        if(refill_counter_inter(rf,w)>= refill_freq(rf,w)){
          
          refill_counter_inter(rf,w) = 0;
          
          if(eq_quantity=="fixed") {
            
            //N_mat(rf,w) =  N_mat(rf,w) + Npu_mat(rf,w);
            N_mat_inter(rf,w) =  N_mat_inter(rf,w)+refill_quantities(rf,w);
            C_mat_inter(rf,w)=0;
            Npu_mat_inter(rf,w)=0;
            
          }
          
        }
        
        else {refill_counter_inter(rf,w)=refill_counter_inter(rf,w)+1;}   
        
        
        
      }
      
    }
    
    if(output== "all") { 
      available_mat_time.slice(time+1) = N_mat;
      cont_mat_time.slice(time+1) = C_mat;
      used_mat_time.slice(time+1) = Npu_mat;
    }
    
  }
  
  //removing last column because because size = time+1
  pop_hosp_status.shed_col(t);
  
  
  if(output=="simple") {

    List L = List::create(_["incidence"] = incidence,
                          _["incidence_inter"] = incidence_inter,
                          _["s_patients"] = s_patients,
                          _["s_patients_inter"] = s_patients_inter,
                          _["cont_new_patients"] = cont_new_patient_vec,
                          _["ward_event"] = ward_event,
                          _["ward_event_inter"] = ward_event_inter,
                          _["count_patients"] = count_patients,
                          _["cont_mat"] = cont_mat,
                          _["used_mat"] = used_mat,
                          _["cont_mat_inter"] = cont_mat_inter,
                          _["used_mat_inter"] = used_mat_inter,
                          _["inf_eq"] = inf_eq,
                          _["inf_eq_inter"] = inf_eq_inter,
                          _["count_patient_ward"] = count_patient_ward,
                          _["count_tests"]= count_tests
  
    );
    
    return L;}
  
  
}
