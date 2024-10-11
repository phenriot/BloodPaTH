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

// Loading all needed functions from these packages
Function random_truncnorm=pkg_truncnorm["rtruncnorm"];
Function random_pert = pkg_mc2d["rpert"];
Function random_lnorm= pkg_stats["rlnorm"];
Function random_norm= pkg_stats["rnorm"];
Function random_unif= pkg_stats["runif"];

//Loading a set_seed function to be used within the model
void set_seed(unsigned int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);  
}

//Data transformation : Returns a list including all associations between procedures and devices to be used within the model's function  
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
List BloodPaTH_model(      float time_step,
                           int t, //simulation time
                           std::string pathogen,
                           int min_e_phase, // min value for the eclipse phase
                           int max_e_phase, // max value for the eclipse phase
                           int nb_patients, // number of patients
                           int nb_wards, // number of wards
                           int nb_adm, // number of admission routes
                           arma::dcolvec adm_prob,// probabilities of admission for each route
                           arma::dmat init_prob, // matrix of probabilities of admission for first ward for both adm routes, nb_adm rows .
                           std::string prev_type, //"ward" level or "admission_route" level
                           arma::drowvec prev_init, // initial prevalence for each admission route
                           arma::dmat WT_matrix, //transition matrix, specific format (merged matrix for all admission routes) but any time-step will do
                           int nb_procedures, //number of procedures
                           StringVector procedure_names, // vector of procedure names
                           int nb_devices, // number of devices
                           StringVector device_names, // vector of devices names
                           arma::imat nb_devices_new, // quantity of device for each type in  each ward
                           arma::imat nb_devices_used, // quantity of device for each type in  each ward
                           arma::imat nb_devices_contaminated, // quantity of device for each type in  each ward
                           arma::imat refill_quantities, // quantity of device to add to the pool of unused device for each type at each refill event
                           arma::fmat refill_freq,// refill frequency for each device type (must be expressed in hours)
                           arma::imat table_procedures_devices, // correspondence table between procedures (col 1) and devices (col 2)
                           arma::drowvec sterilization_prob, // probability of good sterilization for each device
                           arma::dmat PPM_matrix, // probabilities of undergoing procedures (columns) in each ward (row), merged for both admission routes; Any time-step is ok but need to be the same as the transition matrix
                           Rcpp::DataFrame dist_risk, // distribution parameters of the risk for each procedure (dist can be normal, log-normal or pert)
                           std::string intervention = "none",//none, patient-based or ward-based
                           float prob_screening=0, // probability of systematic screening upon admission
                           arma::irowvec at_risk_wards=0, // vector giving the most at risk wards
                           std::string output="simple",//if output == "all" then returns all tables, if "simple" returns nb of new cases each day and new entry in the hospital
                           int id_sim=1
)
{
  
  int real_t = round(t*(1/time_step));
  
  List list_proc_equip = transform_list(nb_procedures,table_procedures_devices);
  
  
  arma::dmat storage_adm_TM;
  arma::dmat storage_adm_PPM;
  arma::drowvec storage_init_prob;
  
  StringVector subset_col_proc = table_procedures_devices[0];
  StringVector subset_col_equip = table_procedures_devices[0];
  
  //Vector of indices to recover each individual transition matrix as the entered matrix is of size ((nb_wards+1)*nb_adm)*(nb_wards+1)
  arma::irowvec index_TM = arma::linspace<arma::irowvec>(0,(nb_wards+1)*nb_adm,nb_adm+1);
  
  //Vector of indices to recover each individual transition matrix as the entered matrix is of size ((nb_wards+1)*nb_adm)*(nb_wards+1)
  arma::irowvec index_PPM = arma::linspace<arma::irowvec>(0,(nb_wards)*nb_adm,nb_adm+1);
  
  
  if(prev_type=="ward" && prev_init.size() != nb_wards) {stop("Error: Initial prevalence type chosen is at a ward-level, should be the same size as the number of wards");}
  if(prev_type=="admission route" && prev_init.size() != nb_adm) {stop("Error: Initial prevalence type chosen is at a admission route-level, should be the same size as the number of admission routes");}
  //if(WT_matrix.n_rows < n_wards) {stop("Error: number of rows in WT_matrix < n_wards");}
  
  //vector for 2= positive, 0=negative
  arma::drowvec pos_neg_vec = {0,2};
  arma::dcolvec cont_neg_vec = {0,1};
  arma::dcolvec failure_success = {0,1};
  arma::drowvec screened_pos_neg_vec = {1,2};
  
  //Matrix for status of patients in the hospital
  arma::dmat pop_hosp_status = arma::dmat(nb_patients,real_t+1,arma::fill::zeros); //No need to t+1 because first index = 0
  arma::dmat pop_hosp_status_inter = arma::dmat(nb_patients,real_t+1,arma::fill::zeros); //No need to t+1 because first index = 0
  
  Rcout << pop_hosp_status.size();
  
  //Matrix for available devices through time
  
  arma::icube available_mat_time(nb_devices,nb_wards,real_t+1);
  arma::icube available_mat_time_inter(nb_devices,nb_wards,real_t+1);
  
  arma::icube cont_mat_time(nb_devices,nb_wards,real_t+1);
  arma::icube cont_mat_time_inter(nb_devices,nb_wards,real_t+1);
  
  arma::icube used_mat_time(nb_devices,nb_wards,real_t+1);
  arma::icube used_mat_time_inter(nb_devices,nb_wards,real_t+1); 
  
  
  // counting nb of patients hospitalized on the periode 
  int count_patients = 0;
  
  // counting nb of tests
  int count_tests = 0;
  
  // counting nb of patients hospitalized on the periode in each ward
  arma::irowvec count_patient_ward(nb_wards,arma::fill::zeros);
  
  // counting nb of contaminated patients hospitalized on the periode in each ward
  arma::irowvec cont_new_patient_ward(nb_wards,arma::fill::zeros);
  
  //Matrix for location of patients in the hospital
  arma::dmat pop_hosp_loc = arma::dmat(nb_patients,real_t+1,arma::fill::zeros); //No need to t+1 because first index = 0
  
  //Matrix of procedures undergone by patients during hosp
  arma::dmat pop_proc = arma::dmat(nb_patients,real_t+1,arma::fill::zeros); //No need to t+1 because first index = 0
  
  //Matrix of number of previously used devicess (initialized at 0 for each procedure in each ward)
  arma::imat Npu_mat = nb_devices_used;
  arma::imat Npu_mat_inter = nb_devices_used;
  
  //Matrix of number of available devices (initialized at nb_devices because no devices has been used at start)
  arma::imat N_mat = nb_devices_new;
  arma::imat N_mat_inter = nb_devices_new;
  
  //Matrix of number of contaminated devicess for each ward and each procedure (C compartment)
  arma::imat C_mat = nb_devices_contaminated;
  arma::imat C_mat_inter = nb_devices_contaminated;
  
  //Matrix storing all departures of patients
  arma::umat patient_departure; //= arma::dmat(2,1,arma::fill::zeros);
  
  //sequence for wards
  arma::drowvec wards_seq = arma::linspace<arma::drowvec>(1,nb_wards,nb_wards);
  
  //sequence for admission places
  arma::drowvec adm_seq = arma::linspace<arma::drowvec>(1,nb_adm,nb_adm);
  
  //sequence for proc
  arma::drowvec proc_seq = arma::linspace<arma::drowvec>(1,nb_procedures,nb_procedures);
  
  //incidence vector
  arma::drowvec incidence(real_t,arma::fill::zeros);
  arma::drowvec incidence_inter(real_t,arma::fill::zeros);
  
  //number of s patients vector
  arma::drowvec s_patients(real_t,arma::fill::zeros);
  arma::drowvec s_patients_inter(real_t,arma::fill::zeros);
  
  //new patients
  arma::drowvec new_patient_vec(real_t,arma::fill::zeros);
  new_patient_vec(0) = nb_patients;
  
  //contaminated new patients
  arma::drowvec cont_new_patient_vec(real_t,arma::fill::zeros);
  arma::drowvec cont_new_patient_vec_inter(t,arma::fill::zeros);
  
  //time counter to assign patient in the E state to the I state
  arma::drowvec state_counter(nb_patients,arma::fill::zeros);// starts at 0
  arma::drowvec state_counter_inter(nb_patients,arma::fill::zeros);// starts at 0
  
  //time counter to assign patient in the E state to the I state
  arma::imat refill_counter(nb_devices,nb_wards,arma::fill::zeros);// starts at 0
  arma::imat refill_counter_inter(nb_devices,nb_wards,arma::fill::zeros);// starts at 0
  
  //store device usage
  arma::imat eq_usage(nb_devices,nb_wards,arma::fill::zeros);// starts at 0
  arma::imat eq_usage_inter(nb_devices,nb_wards,arma::fill::zeros);// starts at 0
  
  //time counter to assign patient in the E state to the I state
  arma::drowvec admission_route(nb_patients,arma::fill::zeros);// starts at 0
  //arma::drowvec admission_route_inter(nb_patients,arma::fill::zeros);// starts at 0
  
  //screened vector upon admission
  arma::drowvec screened(nb_patients,arma::fill::zeros);// starts at 0
  
  arma::dcolvec screen_prob_vec = {1-prob_screening,prob_screening};
  
  //Matrix to know in patient has been screened for HCV or HBV or HIV , 0 = not screened, 1 = screened and - , 2 = screened and +
  arma::imat ward_screening = arma::imat(nb_patients,nb_wards,arma::fill::zeros); //No need to t+1 because first index = 0
  
  //vector to store nb of infection events
  arma::drowvec ward_event(nb_wards,arma::fill::zeros);// starts at 0
  arma::drowvec ward_event_inter(nb_wards,arma::fill::zeros);// starts at 0
  
  //matrix to store nb of suceptible patients in each ward
  arma::dmat ward_s = arma::dmat(nb_wards,real_t+1,arma::fill::zeros);
  arma::dmat ward_s_inter = arma::dmat(nb_wards,real_t+1,arma::fill::zeros);
  
  //vector to store nb of infection events
  arma::drowvec cont_mat(nb_devices,arma::fill::zeros);// starts at 0
  arma::drowvec cont_mat_inter(nb_devices,arma::fill::zeros);// starts at 0
  
  //vector to store nb of infection events
  arma::drowvec used_mat(nb_devices,arma::fill::zeros);// starts at 0
  arma::drowvec used_mat_inter(nb_devices,arma::fill::zeros);// starts at 0
  
  //vector to store nb of proc by proc
  arma::drowvec vec_nb_proc_inter(nb_procedures,arma::fill::zeros);// starts at 0
  arma::drowvec vec_nb_proc(nb_procedures,arma::fill::zeros);// starts at 0
  
  //vector to store nb of infection per device type
  arma::drowvec inf_eq(nb_devices,arma::fill::zeros);// starts at 0
  arma::drowvec inf_eq_inter(nb_devices,arma::fill::zeros);// starts at 0
  
  //Matrix to know if a patient has already been into a given ward
  arma::imat mat_patient_ward = arma::imat(nb_patients,nb_wards,arma::fill::zeros); //No need to t+1 because first index = 0
  
  arma::irowvec init_patient_wards(nb_wards,arma::fill::zeros);
  
  arma::irowvec fill_screened_s = std::vector<int>(nb_wards,1);
  
  arma::irowvec fill_screened_i = std::vector<int>(nb_wards,2);
  
  
  //
  arma::drowvec e_phase = Rcpp::runif( nb_patients, round(min_e_phase*(1/time_step)), round(max_e_phase*(1/time_step)));
  e_phase = round(e_phase);
  
  arma::drowvec patient_id =  arma::linspace<arma::drowvec>(0,nb_patients-1,nb_patients);
  set_seed(1*id_sim);
  arma::drowvec rdm_patient  = sample(patient_id,nb_patients,false);
  
  // For time t=0 initial status are added to pop_hosp
  
  for(int i = 0; i<nb_patients;i++) {
    
    
    //set_seed(id_sim); //on reset les seed avant chaque sample pour les differents patients pour s'assurer de la repro de certains processus
    admission_route(i) = sample(adm_seq,1,false,adm_prob)(0);
    //Rcout<<admission_ward(i);
    storage_init_prob= init_prob.row(admission_route(i)-1);
    //Rcout<<wards_seq;
    
    pop_hosp_loc(i,0) = sample(wards_seq,1,false,storage_init_prob.t())(0);
    //Rcout << "ici";
    
    if(prev_type == "admission route") {
      
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
  
  for(unsigned int time = 0; time<real_t;time++) {
    
    rdm_patient  = sample(patient_id,nb_patients,false);
    
    for(unsigned int j = 0; j<nb_patients;j++) {
      
      unsigned int p = rdm_patient(j);
      
      int current_ward = pop_hosp_loc(p,time);

      // giving submatrix : rows in span  index_TM(0) --> index_TM(1)-1
      storage_adm_TM = WT_matrix(arma::span(index_TM(admission_route(p)-1),index_TM(admission_route(p))-1),arma::span(0,nb_wards));
      storage_adm_TM =  storage_adm_TM.t();
      arma::dcolvec prob_ward =storage_adm_TM.col(current_ward-1);//minus 1 because first index = 0
      
      storage_adm_PPM = PPM_matrix(arma::span(index_PPM(admission_route(p)-1),index_PPM(admission_route(p))-1),arma::span(0,(nb_procedures-1)));
      storage_adm_PPM =storage_adm_PPM.t();
      
      pop_proc(p,time) = sample(proc_seq,1,false, storage_adm_PPM.col(pop_hosp_loc(p,time)-1))(0); // random draw of proc in prob linked to selected ward
      
      
      current_ward = pop_hosp_loc(p,time);
      
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
      
      if(time>0) {
        if((intervention == "ward-based") && 
           screened(p) == 0 && 
           std::find(at_risk_wards.begin(), at_risk_wards.end(), current_ward) != at_risk_wards.end() &&
           pop_hosp_loc(p,time-1) != pop_hosp_loc(p,time)) {
          
          int screen = sample(failure_success,1,false,screen_prob_vec)(0);
          
          if (screen==1 && pop_hosp_status_inter(p,time-1)==0) {screened(p) = 1;  count_tests = count_tests +1 ;}
          if (screen==1 && pop_hosp_status_inter(p,time-1)==2) {screened(p) = 2;  count_tests = count_tests +1 ;}
          
        } 
        
      }
      
      
      //////////////////////////
      // NO INTERVENTION PART //
      //////////////////////////
      
      //for infected patients
      //current_proc != nb_proc because nb_procth procedure = no procedure
      if((viremic_status == 2 or viremic_status == 3) && current_proc != nb_procedures && current_ward != state_out) {
        
        
        //list of involved devices in given procedure
        arma::irowvec list_equip = list_proc_equip(current_proc-1);
        
        //loop over involved devices
        for(int eq = 0; eq< list_equip.n_cols; eq++ ) {
          
          int current_device = list_equip(eq); 
          
          eq_usage(current_device-1,current_ward-1) = eq_usage(current_device-1,current_ward-1) + 1;
          
          float n_available = N_mat(current_device-1,current_ward-1); //nb of total available devices for a given proc and ward
          
          float n_pu = Npu_mat(current_device-1,current_ward-1); //nb of total previously used for a given proc and ward
          
          
          //random draw for reused devices
          int reuse = 0;
          
          
          float init_eq = n_available + n_pu;
          
          //if no more available devices
          if(n_available < 1) {reuse = 1;}
          
          //if no previously used devices
          if(n_pu == 0) {reuse = 0;}
          
          //probability of reused devices being contaminated
          
          //if reusing draw to know if reused mat is contaminated
          if(reuse == 1) {
            
            used_mat(current_device-1) = used_mat(current_device-1)+1;
            
            float n_c = C_mat(current_device-1,current_ward-1);
            double prob_Cpw = n_c / n_pu; // prob of drawing a contaminated device
            
            if(prob_Cpw > 1) {Rcout<< prob_Cpw; prob_Cpw=1;}
            
            arma::dcolvec vec_prob_Cpw = {1-prob_Cpw,prob_Cpw};
            
            set_seed(round((((time+1)+(p+1))/(p+1))*id_sim*1));
            int draw_mat = sample(failure_success,1,false,vec_prob_Cpw)(0);
            
            //if devices not previously contaminated, it becomes contaminated with a probability equal to sterilization efficiency for the given device
            if(draw_mat==0){
              
              // probability of device being well sterilized
              double prob_ster = sterilization_prob(current_device-1);
              arma::dcolvec vec_prob_ster = {1-prob_ster,prob_ster};
              
              set_seed(round((((time+1)+(p+1))/(p+1))*id_sim*2));
              int sterilized = sample(failure_success,1,false,vec_prob_ster)(0);
              
              //if patient is detected as positive
              // if sterilization fails
              if(sterilized==0) {C_mat(current_device-1,current_ward-1) = C_mat(current_device-1,current_ward-1)+1;}
              
            }
            
            //if devices previously contaminated, it can be sterilized
            if(draw_mat==1){
              
              cont_mat(current_device-1) = cont_mat(current_device-1)+1;
              
              // probability of device being well sterilized
              double prob_ster = sterilization_prob(current_device-1);
              arma::dcolvec vec_prob_ster = {1-prob_ster,prob_ster};
              
              set_seed(round((((time+1)+(p+1))/(p+1))*id_sim*3));
              int sterilized = sample(failure_success,1,false,vec_prob_ster)(0);
              
              // if sterilization is effective
              if(sterilized==1) {C_mat(current_device-1,current_ward-1) = C_mat(current_device-1,current_ward-1)-1;}
              
            }
            
          }
          
          //if not reusing devices then it gets contaminated depending on the good sterilization prob
          if(reuse == 0) {
            
            //adding 1 to matrix of previously used devices
            Npu_mat(current_device-1,current_ward-1) = Npu_mat(current_device-1,current_ward-1) + 1;
            
            //removing 1 to matrix of available devices
            N_mat(current_device-1,current_ward-1) = N_mat(current_device-1,current_ward-1)-1;
            
            double prob_ster = sterilization_prob(current_device-1);
            arma::dcolvec vec_prob_ster = {1-prob_ster,prob_ster};
            
            set_seed(round((((time+1)+(p+1))/(p+1))*id_sim*4));
            int sterilized = sample(failure_success,1,false,vec_prob_ster)(0);
            
            
            //devices becomes contaminated
            if(sterilized==0) {C_mat(current_device-1,current_ward-1) = C_mat(current_device-1,current_ward-1)+1;}
            
          }
        }
      }
      
      //patients whom are infected but not yet infectious
      if (viremic_status == 1 && current_proc != nb_procedures && current_ward != state_out) {
        
        //list of involved devices in given procedure
        arma::irowvec list_equip = list_proc_equip(current_proc-1);
        
        //loop over involved devices
        for(int eq = 0; eq< list_equip.n_cols; eq++ ) {
          
          int current_device = list_equip(eq);
          
          eq_usage(current_device-1,current_ward-1) = eq_usage(current_device-1,current_ward-1) + 1;
          
          float n_available = N_mat(current_device-1,current_ward-1); //nb of total available devices for a given proc and ward
          
          float n_pu = Npu_mat(current_device-1,current_ward-1); //nb of total previously used for a given proc and ward
          
          
          //random draw for reused devices
          int reuse = 0;
          
          float init_eq = n_available + n_pu;
          
          //if no more available devices
          if(n_available < 1) {reuse = 1;}
          
          
          //if no previously used devices
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
            double prob_Cpw = n_c / n_pu; // prob of drawing a contaminated device
            
            if(prob_Cpw > 1) {Rcout<< prob_Cpw; prob_Cpw=1;}
            
            arma::dcolvec vec_prob_Cpw = {1-prob_Cpw,prob_Cpw};
            
            set_seed(round((((time+1)+(p+1))/(p+1))*id_sim*5));
            int draw_mat = sample(failure_success,1,false,vec_prob_Cpw)(0);
            
            //if devices contaminated
            if(draw_mat==1){
              
              cont_mat(current_device-1) = cont_mat(current_device-1)+1;
              
              // probability of device being well sterilized
              double prob_ster = sterilization_prob(current_device-1);
              arma::dcolvec vec_prob_ster = {1-prob_ster,prob_ster};
              
              set_seed(round((((time+1)+(p+1))/(p+1))*id_sim*6));
              int sterilized = sample(failure_success,1,false,vec_prob_ster)(0);
              
              // if sterilization is well performed
              if(sterilized==1) {C_mat(current_device-1,current_ward-1) = C_mat(current_device-1,current_ward-1)-1;}
              
            }
          }
        }
      }
      
      
      //susceptible patients
      if (viremic_status == 0 && current_proc != nb_procedures && current_ward != state_out) {
        
        //list of involved devices in given procedure
        arma::irowvec list_equip = list_proc_equip(current_proc-1);
        
        // Vector of the prob of device being contaminated (used if multiple devices for one procedure)
        arma::drowvec prob_eq_vec(list_equip.n_cols);
        
        //loop over involved devices
        for(int eq = 0; eq< list_equip.n_cols; eq++ ) {
          
          int current_device = list_equip(eq); 
          
          eq_usage(current_device-1,current_ward-1) = eq_usage(current_device-1,current_ward-1) + 1;
          
          float n_available = N_mat(current_device-1,current_ward-1); //nb of total available devices for a given proc and ward
          
          float n_pu = Npu_mat(current_device-1,current_ward-1); //nb of total previously used for a given proc and ward
          
          //random draw for reused devices
          int reuse = 0;
          
          //if no more available devices
          float init_eq = n_available + n_pu;
          
          if(n_available < 1) {reuse = 1;}
          
          // if(n_available > 1 && n_available/init_eq<threshold_reuse(current_device-1)) {
          //   
          //   double prob_reuse= 1-(n_available/init_eq);
          //   
          //   arma::dcolvec vec_prob_reuse = {1-prob_reuse,prob_reuse};
          //   
          //   set_seed(round((((time+1)+(p+1))/(p+1))*id_sim*7));
          //   reuse = sample(failure_success,1,false,vec_prob_reuse)(0);}
          
          //if no previously used devices
          if(n_pu == 0) {reuse = 0;}
          
          //probability of reused devices being contaminated
          
          //if reusing draw to know if reused mat is contaminated
          if(reuse == 1) {
            
            used_mat(current_device-1) = used_mat(current_device-1)+1;
            
            float n_c = C_mat(current_device-1,current_ward-1);
            double prob_Cpw = n_c / n_pu; // prob of drawing a contaminated device
            
            if(prob_Cpw > 1) {Rcout<< prob_Cpw; prob_Cpw=1;}
            
            arma::dcolvec vec_prob_Cpw = {1-prob_Cpw,prob_Cpw};
            
            set_seed(round((((time+1)+(p+1))/(p+1))*id_sim*8));
            int draw_mat = sample(failure_success,1,false,vec_prob_Cpw)(0);
            
            prob_eq_vec(eq) = draw_mat;
            
            //if devices contaminated
            if(draw_mat==1){
              
              cont_mat(current_device-1) = cont_mat(current_device-1)+1;
              
              // probability of device being well sterilized
              double prob_ster = sterilization_prob(current_device-1);
              arma::dcolvec vec_prob_ster = {1-prob_ster,prob_ster};
              
              set_seed(round((((time+1)+(p+1))/(p+1))*id_sim*9));
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
        
        //if a least 1 eq is contaminated
        if(sum(prob_eq_vec) >= 1) {
          
          
          StringVector vec_dist = dist_risk["dist"];
          
          if(vec_dist(current_proc-1)=="lnorm") {
            
            arma::fvec lmean_vec = dist_risk[0];
            
            arma::fvec lsd_vec = dist_risk[1];
            
            float lmean_risk = lmean_vec(current_proc-1); //parameter 1 for lognormal dist
            
            float lsd_risk = lsd_vec(current_proc-1); //parameter 2 for lognormal dist
            
            set_seed(round((((time+1)+(p+1))/(p+1))*id_sim*10));
            float risk_p = Rcpp::as<float> (random_lnorm(_["n"]=1,_["meanlog"]=lmean_risk,_["sdlog"]=lsd_risk)); //risk of getting HCV contaminated if devices is contaminated
            
            float prob_E =risk_p ; 
            
            if(pathogen =="HBV") {
              
              set_seed(round((((time+1)+(p+1))/(p+1))*id_sim*11));
              float HBV_mult = Rcpp::as<float> (random_unif(_["n"]=1,_["min"]=4,_["max"]=20)); //risk mult for HBV
              
              prob_E = risk_p*HBV_mult;
              
            }
            
            
            //Probability of going from the S to the E state
            
            if (prob_E<0) {prob_E=0;}
            if (prob_E>1) {prob_E=1;}
            
            arma::dcolvec prob_inf_vector = {1-prob_E,prob_E};
            
            set_seed(round((((time+1)+(p+1))/(p+1))*id_sim*12));
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
            
            set_seed(round((((time+1)+(p+1))/(p+1))*id_sim*13));
            float risk_p = Rcpp::as<float> (random_norm(_["n"]=1,_["mean"]=mean_risk,_["sd"]=sd_risk)); //risk of getting HCV contaminated if devices is contaminated
            
            float prob_E = risk_p;
            
            if(pathogen=="HBV") {
              
              set_seed(round((((time+1)+(p+1))/(p+1))*id_sim*14));
              float HBV_mult = Rcpp::as<float> (random_unif(_["n"]=1,_["min"]=4,_["max"]=20)); //risk mult for HBV
              
              prob_E = risk_p*HBV_mult;
              
            }
            
            if (prob_E<0) {prob_E=0;}
            if (prob_E>1) {prob_E=1;}
            
            arma::dcolvec prob_inf_vector = {1-prob_E,prob_E};
            
            set_seed(round((((time+1)+(p+1))/(p+1))*id_sim*15));
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
            
            set_seed(round((((time+1)+(p+1))/(p+1))*id_sim*16));
            float risk_p = Rcpp::as<float> (random_pert(_["n"]=1,_["min"]=min_risk,_["mode"]=mode_risk,_["max"]=max_risk))/100; //risk of getting HCV contaminated if devices is contaminated
            
            float prob_E = risk_p;
            
            if(pathogen =="HBV") {
              
              set_seed(round((((time+1)+(p+1))/(p+1))*id_sim*17));
              float HBV_mult = Rcpp::as<float> (random_unif(_["n"]=1,_["min"]=4,_["max"]=20)); //risk mult for HBV
              
              prob_E = risk_p*HBV_mult;
              
            }
            
            if (prob_E<0) {prob_E=0;}
            if (prob_E>1) {prob_E=1;}
            
            arma::dcolvec prob_inf_vector = {1-prob_E,prob_E};
            
            set_seed(round((((time+1)+(p+1))/(p+1))*id_sim*18));
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
      
      //for infected patients
      //current_proc != nb_proc because nb_procth procedure = no procedure
      if((viremic_status_inter == 2 or viremic_status_inter == 3) && current_proc != nb_procedures && current_ward != state_out) {
        
        
        //list of involved devices in given procedure
        arma::irowvec list_equip = list_proc_equip(current_proc-1);
        
        //loop over involved devices
        for(int eq = 0; eq< list_equip.n_cols; eq++ ) {
          
          int current_device = list_equip(eq); 
          
          eq_usage_inter(current_device-1,current_ward-1) = eq_usage_inter(current_device-1,current_ward-1) + 1;
          
          float n_available = N_mat_inter(current_device-1,current_ward-1); //nb of total available devices for a given proc and ward
          
          float n_pu = Npu_mat_inter(current_device-1,current_ward-1); //nb of total previously used for a given proc and ward
          
          //random draw for reused devices
          int reuse = 0;
          
          float init_eq = n_available + n_pu;
          
          //if no more available devices
          if(n_available < 1) {reuse = 1;}
          
          //if no previously used devices
          if(n_pu == 0) {reuse = 0;}
          
          //probabiliused_matty of reused devices being contaminated
          
          //if reusing draw to know if reused mat is contaminated
          if(reuse == 1) {
            
            used_mat_inter(current_device-1) = used_mat_inter(current_device-1)+1;
            
            float n_c = C_mat_inter(current_device-1,current_ward-1);
            double prob_Cpw = n_c / n_pu; // prob of drawing a contaminated device
            
            arma::dcolvec vec_prob_Cpw = {1-prob_Cpw,prob_Cpw};
            
            set_seed(round((((time+1)+(p+1))/(p+1))*id_sim*1));
            int draw_mat = sample(failure_success,1,false,vec_prob_Cpw)(0);
            
            //if devices not previously contaminated, it becomes contaminated with a probability equal to sterilization efficiency for the given device
            if(draw_mat==0){
              
              // probability of device being well sterilized
              double prob_ster = sterilization_prob(current_device-1);
              arma::dcolvec vec_prob_ster = {1-prob_ster,prob_ster};
              
              set_seed(round((((time+1)+(p+1))/(p+1))*id_sim*2));
              int sterilized = sample(failure_success,1,false,vec_prob_ster)(0);
              
              //if patient is detected as positive
              if( (intervention == "patient-based" or intervention == "ward-based") && screened(p) == 2) {sterilized = 1;}
              
              // if sterilization fails
              if(sterilized==0) {C_mat_inter(current_device-1,current_ward-1) = C_mat_inter(current_device-1,current_ward-1)+1;}
              
            }
            
            //if devices previously contaminated, it can be sterilized
            if(draw_mat==1){
              
              cont_mat_inter(current_device-1) = cont_mat_inter(current_device-1)+1;
              
              // probability of device being well sterilized
              double prob_ster = sterilization_prob(current_device-1);
              arma::dcolvec vec_prob_ster = {1-prob_ster,prob_ster};
              
              set_seed(round((((time+1)+(p+1))/(p+1))*id_sim*3));
              int sterilized = sample(failure_success,1,false,vec_prob_ster)(0);
              
              if( (intervention == "patient-based" or intervention == "ward-based") && screened(p) == 2) {sterilized = 1;}
              
              // if sterilization is effective
              if(sterilized==1) {C_mat_inter(current_device-1,current_ward-1) = C_mat_inter(current_device-1,current_ward-1)-1;}
              
            }
          }
          
          //if not reusing devices then it gets contaminated depending on the good sterilization prob
          if(reuse == 0) {
            
            //adding 1 to matrix of previously used devices
            Npu_mat_inter(current_device-1,current_ward-1) = Npu_mat_inter(current_device-1,current_ward-1) + 1;
            
            //removing 1 to matrix of available devices
            N_mat_inter(current_device-1,current_ward-1) = N_mat_inter(current_device-1,current_ward-1)-1;
            
            double prob_ster = sterilization_prob(current_device-1);
            arma::dcolvec vec_prob_ster = {1-prob_ster,prob_ster};
            
            set_seed(round((((time+1)+(p+1))/(p+1))*id_sim*4));
            int sterilized = sample(failure_success,1,false,vec_prob_ster)(0);
            
            if( (intervention == "patient-based" or intervention == "ward-based") && screened(p) == 2) {sterilized = 1;}
            
            //devices becomes contaminated
            if(sterilized==0) {C_mat_inter(current_device-1,current_ward-1) = C_mat_inter(current_device-1,current_ward-1)+1;}
            
          }
        }
      }
      
      //patients that are infected but not infectious
      if (viremic_status_inter == 1 && current_proc != nb_procedures && current_ward != state_out) {
        
        //list of involved devices in given procedure
        arma::irowvec list_equip = list_proc_equip(current_proc-1);
        
        //loop over involved devices
        for(int eq = 0; eq< list_equip.n_cols; eq++ ) {
          
          int current_device = list_equip(eq);
          
          eq_usage_inter(current_device-1,current_ward-1) = eq_usage_inter(current_device-1,current_ward-1) + 1;
          
          float n_available = N_mat_inter(current_device-1,current_ward-1); //nb of total available devices for a given proc and ward
          
          float n_pu = Npu_mat_inter(current_device-1,current_ward-1); //nb of total previously used for a given proc and ward
          
          
          //random draw for reused devices
          int reuse = 0;
          
          float init_eq = n_available + n_pu;
          //Rcout << init_eq;
          
          //if no more available devices
          if(n_available < 1) {reuse = 1;}
          
          //if no previously used devices
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
            double prob_Cpw = n_c / n_pu; // prob of drawing a contaminated device
            
            if(prob_Cpw > 1) {Rcout<< prob_Cpw; prob_Cpw=1;}
            
            arma::dcolvec vec_prob_Cpw = {1-prob_Cpw,prob_Cpw};
            
            set_seed(round((((time+1)+(p+1))/(p+1))*id_sim*5));
            int draw_mat = sample(failure_success,1,false,vec_prob_Cpw)(0);
            
            //if devices contaminated
            if(draw_mat==1){
              
              cont_mat_inter(current_device-1) = cont_mat_inter(current_device-1)+1;
              
              // probability of device being well sterilized
              double prob_ster = sterilization_prob(current_device-1);
              arma::dcolvec vec_prob_ster = {1-prob_ster,prob_ster};
              
              set_seed(round((((time+1)+(p+1))/(p+1))*id_sim*6));
              int sterilized = sample(failure_success,1,false,vec_prob_ster)(0);
              
              // if sterilization is well performed
              if(sterilized==1) {C_mat_inter(current_device-1,current_ward-1) = C_mat_inter(current_device-1,current_ward-1)-1;}
              
            }
          }
        }
      }
      
      //susceptible patients
      if (viremic_status_inter == 0 && current_proc != nb_procedures && current_ward != state_out) {
        
        //list of involved devices in given procedure
        arma::irowvec list_equip = list_proc_equip(current_proc-1);
        
        // Vector of the prob of device being contaminated (used if multiple devices for one procedure)
        arma::drowvec prob_eq_vec(list_equip.n_cols);
        
        //loop over involved devices
        for(int eq = 0; eq< list_equip.n_cols; eq++ ) {
          
          int current_device = list_equip(eq); 
          
          eq_usage_inter(current_device-1,current_ward-1) = eq_usage_inter(current_device-1,current_ward-1) + 1;
          
          float n_available = N_mat_inter(current_device-1,current_ward-1); //nb of total available devices for a given proc and ward
          
          float n_pu = Npu_mat_inter(current_device-1,current_ward-1); //nb of total previously used for a given proc and ward
          
          //random draw for reused devices
          int reuse = 0;
          
          //if no more available devices
          float init_eq = n_available + n_pu;
          
          if(n_available < 1) {reuse = 1;}
          
          // if(n_available > 1 && n_available/init_eq<threshold_reuse(current_device-1)) {
          //   
          //   double prob_reuse= 1-(n_available/init_eq);
          //   
          //   arma::dcolvec vec_prob_reuse = {1-prob_reuse,prob_reuse};
          //   
          //   //Rcout << vec_prob_reuse;
          //   
          //   set_seed(round((((time+1)+(p+1))/(p+1))*id_sim*7));
          //   reuse = sample(failure_success,1,false,vec_prob_reuse)(0);}
          
          //if no previously used devices
          if(n_pu == 0) {reuse = 0;}
          
          
          //probability of reused devices being contaminated
          
          //if reusing draw to know if reused mat is contaminated
          if(reuse == 1) {
            
            used_mat_inter(current_device-1) = used_mat_inter(current_device-1)+1;
            
            float n_c = C_mat_inter(current_device-1,current_ward-1);
            double prob_Cpw = n_c / n_pu; // prob of drawing a contaminated device
            
            arma::dcolvec vec_prob_Cpw = {1-prob_Cpw,prob_Cpw};
            
            set_seed(round((((time+1)+(p+1))/(p+1))*id_sim*8));
            int draw_mat = sample(failure_success,1,false,vec_prob_Cpw)(0);
            
            prob_eq_vec(eq) = draw_mat;
            
            //if devices contaminated
            if(draw_mat==1){
              
              cont_mat_inter(current_device-1) = cont_mat_inter(current_device-1)+1;
              
              // probability of device being well sterilized
              double prob_ster = sterilization_prob(current_device-1);
              arma::dcolvec vec_prob_ster = {1-prob_ster,prob_ster};
              
              set_seed(round((((time+1)+(p+1))/(p+1))*id_sim*9));
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
        
        //if a least 1 eq is contaminated
        if(sum(prob_eq_vec) >= 1) {
          
          StringVector vec_dist = dist_risk["dist"];
          
          if(vec_dist(current_proc-1)=="lnorm") {
            
            arma::fvec lmean_vec = dist_risk[0];
            
            arma::fvec lsd_vec = dist_risk[1];
            
            float lmean_risk = lmean_vec(current_proc-1); //parameter 1 for lognormal dist
            
            float lsd_risk = lsd_vec(current_proc-1); //parameter 2 for lognormal dist
            
            set_seed(round((((time+1)+(p+1))/(p+1))*id_sim*10));
            float risk_p = Rcpp::as<float> (random_lnorm(_["n"]=1,_["meanlog"]=lmean_risk,_["sdlog"]=lsd_risk)); //risk of getting HCV contaminated if devices is contaminated
            
            float prob_E =risk_p ; 
            
            if(pathogen =="HBV") {
              
              set_seed(round((((time+1)+(p+1))/(p+1))*id_sim*11));
              float HBV_mult = Rcpp::as<float> (random_unif(_["n"]=1,_["min"]=4,_["max"]=20)); //risk mult for HBV
              
              prob_E = risk_p*HBV_mult;
              
            }
            
            //Probability of going from the S to the E state
            
            if (prob_E<0) {prob_E=0;}
            if (prob_E>1) {prob_E=1;}
            
            arma::dcolvec prob_inf_vector = {1-prob_E,prob_E};
            
            set_seed(round((((time+1)+(p+1))/(p+1))*id_sim*12));
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
            
            set_seed(round((((time+1)+(p+1))/(p+1))*id_sim*13));
            float risk_p = Rcpp::as<float> (random_norm(_["n"]=1,_["mean"]=mean_risk,_["sd"]=sd_risk)); //risk of getting HCV contaminated if devices is contaminated
            
            float prob_E = risk_p;
            
            if(pathogen=="HBV") {
              
              set_seed(round((((time+1)+(p+1))/(p+1))*id_sim*14));
              float HBV_mult = Rcpp::as<float> (random_unif(_["n"]=1,_["min"]=4,_["max"]=20)); //risk mult for HBV
              
              prob_E = risk_p*HBV_mult;
              
            }
            
            if (prob_E<0) {prob_E=0;}
            if (prob_E>1) {prob_E=1;}
            
            arma::dcolvec prob_inf_vector = {1-prob_E,prob_E};
            
            set_seed(round((((time+1)+(p+1))/(p+1))*id_sim*15));
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
            
            set_seed(round((((time+1)+(p+1))/(p+1))*id_sim*16));
            float risk_p = Rcpp::as<float> (random_pert(_["n"]=1,_["min"]=min_risk,_["mode"]=mode_risk,_["max"]=max_risk))/100; //risk of getting HCV contaminated if devices is contaminated
            
            float prob_E = risk_p;
            
            if(pathogen =="HBV") {
              
              set_seed(round((((time+1)+(p+1))/(p+1))*id_sim*17));
              float HBV_mult = Rcpp::as<float> (random_unif(_["n"]=1,_["min"]=4,_["max"]=20)); //risk mult for HBV
              
              prob_E = risk_p*HBV_mult;
              
            }
            
            if (prob_E<0) {prob_E=0;}
            if (prob_E>1) {prob_E=1;}
            
            arma::dcolvec prob_inf_vector = {1-prob_E,prob_E};
            
            
            set_seed(round((((time+1)+(p+1))/(p+1))*id_sim*18));
            int draw_infection = sample(cont_neg_vec,1,false,prob_inf_vector)(0);
            
            if (draw_infection == 1) {pop_hosp_status_inter(p,time) = 1; incidence_inter(time)=incidence_inter(time)+1; ward_event_inter(current_ward-1)= ward_event_inter(current_ward-1)+1;
            
            for(int ei = 0; ei < list_equip.n_cols; ei++) {if(prob_eq_vec(ei)==1) {inf_eq_inter(list_equip(ei)-1) = inf_eq_inter(list_equip(ei)-1)+1;} } 
            
            }
          }
        }
      }
      
      //adding 1 to status counter
      if(pop_hosp_status(p,time) == 1) {state_counter(p) = state_counter(p) + 1;}
      if(pop_hosp_status_inter(p,time) == 1) {state_counter_inter(p) = state_counter_inter(p) + 1;}
      
      
      //if counter exceed the time of the eclipse phase then patient goes from E to I
      if(state_counter(p) >= e_phase(p) && pop_hosp_status(p,time) == 1 ) {pop_hosp_status(p,time) = 3; }
      if(state_counter_inter(p) >= e_phase(p) && pop_hosp_status_inter(p,time) == 1 ) {pop_hosp_status_inter(p,time) = 3; }
      
      pop_hosp_loc(p,time+1) = sample(wards_out_seq,1,false,prob_ward)(0); //first ward for new patient
      
      
      pop_hosp_status(p,time+1) = pop_hosp_status(p,time); //patient has same status than time before
      pop_hosp_status_inter(p,time+1) = pop_hosp_status_inter(p,time); //patient has same status than time before
      
      
      
      if(pop_hosp_status(p,time)==0) {s_patients(time)=s_patients(time)+1; ward_s(current_ward-1,time)=ward_s(current_ward-1,time)+1;}
      if(pop_hosp_status_inter(p,time)==0) {s_patients_inter(time)=s_patients_inter(time)+1; ward_s_inter(current_ward-1,time)=ward_s_inter(current_ward-1,time)+1;}
      
      // if a patient leaves the hospital, another one is entering
      if( pop_hosp_loc(p,time+1) == state_out) {
        
        //initialize to zero so that history of passing through wards is clear for new patient
        mat_patient_ward.row(p) = init_patient_wards;
        
        ward_screening.row(p) = init_patient_wards;
        
        count_patients = count_patients +1;
        
        new_patient_vec(time) = new_patient_vec(time)+1;
        
        arma::ucolvec new_departure = {time+1,p};
        patient_departure = join_horiz(patient_departure ,new_departure) ;
        
        admission_route(p) = sample(adm_seq,1,false,adm_prob)(0);
        storage_init_prob= init_prob.row(admission_route(p)-1);
        pop_hosp_loc(p,time+1) = sample(wards_seq,1,false,storage_init_prob.t())(0);//if patient leaves the hospital, another patient instantly replaces him
        
        state_counter(p) = 0; //if patient leaves the hospital, counter of infection for patient p goes back to 0
        state_counter_inter(p) = 0;
        
        if(prev_type == "ward") {
          
          arma::dcolvec pos_neg_prob = {1-prev_init((pop_hosp_loc(p,time+1)-1)),prev_init((pop_hosp_loc(p,time+1)-1))};
          
          //test if patient is infected depending on adm route
          pop_hosp_status(p,time+1) = sample(pos_neg_vec,1,false,pos_neg_prob)(0);
          pop_hosp_status_inter(p,time+1) = pop_hosp_status(p,time+1);
          
        }
        
        if(prev_type == "admission route") {
          
          arma::dcolvec pos_neg_prob = {1-prev_init((admission_route(p)-1)),prev_init((admission_route(p)-1))};
          
          //test if patient is infected depending on adm route
          pop_hosp_status(p,time+1) = sample(pos_neg_vec,1,false,pos_neg_prob)(0);
          pop_hosp_status_inter(p,time+1) = pop_hosp_status(p,time+1);
          
        }
        
        e_phase(p) = round(Rcpp::runif(1,round(min_e_phase*(1/time_step)),round(max_e_phase*(1/time_step)))(0));
        
        
        screened(p)=0;
        
        if(intervention == "ward-based" && std::find(at_risk_wards.begin(), at_risk_wards.end(), (pop_hosp_loc(p,time+1))) != at_risk_wards.end()) {
          
          int screen = sample(failure_success,1,false,screen_prob_vec)(0);
          
          if (screen==1 && pop_hosp_status_inter(p,time+1)==0) {screened(p)=1; count_tests = count_tests +1 ;}
          if (screen==1 && pop_hosp_status_inter(p,time+1)==2) {screened(p)=2; count_tests = count_tests +1 ;}
          
        }
        
        if(intervention =="patient-based") {
          
          int screen_new = sample(failure_success,1,false,screen_prob_vec)(0);
          
          if (screen_new==1 && pop_hosp_status_inter(p,time+1)==0) {screened(p)=1;  count_tests = count_tests +1 ;}
          
          if (screen_new==1 && pop_hosp_status_inter(p,time+1)==2) {screened(p)=2;  count_tests = count_tests +1 ;}
          
        }
      }
    }
    
    
    // check if there s a need to refill device or throw away device
    for(unsigned int rf=0;rf<nb_devices;rf++) {
      
      for(unsigned int w = 0; w<nb_wards; w++) {
        
        // if(eq_quantity=="variable") {
        //   if(refill_counter(rf,w)>= round(refill_freq(rf,w)*(1/time_step))) {N_mat(rf,w)=refill_quantities(rf,w); refill_counter(rf,w) = 0;}
        //   else {refill_counter(rf,w)=refill_counter(rf,w)+1;}
        //   
        //   if(refill_counter_inter(rf,w)>= round(refill_freq(rf,w)*(1/time_step))) {N_mat_inter(rf,w)=refill_quantities(rf,w); refill_counter_inter(rf,w) = 0;}
        //   else {refill_counter_inter(rf,w)=refill_counter_inter(rf,w)+1;}
        //   
        // }
        // 
        if(refill_counter(rf,w)>= round(refill_freq(rf,w)*(1/time_step))){
          
          refill_counter(rf,w) = 0;
          
          // if(eq_quantity=="fixed") {
          
          N_mat(rf,w) =  N_mat(rf,w)+refill_quantities(rf,w);
          C_mat(rf,w)=0;
          Npu_mat(rf,w)=0;
          
          // }
        }
        
        else {refill_counter(rf,w)=refill_counter(rf,w)+1;}
        
        if(refill_counter_inter(rf,w)>= round(refill_freq(rf,w)*(1/time_step))){
          
          refill_counter_inter(rf,w) = 0;
          
          // if(eq_quantity=="fixed") {
          
          N_mat_inter(rf,w) =  N_mat_inter(rf,w)+refill_quantities(rf,w);
          C_mat_inter(rf,w)=0;
          Npu_mat_inter(rf,w)=0;
          
          // }
        }
        
        else {refill_counter_inter(rf,w)=refill_counter_inter(rf,w)+1;}   
        
      }
    }
    
    // if(output== "all") { 
    //   available_mat_time.slice(time+1) = N_mat;
    //   cont_mat_time.slice(time+1) = C_mat;
    //   used_mat_time.slice(time+1) = Npu_mat;
    //   
    // }
  }
  
  //removing last column because because size = time+1
  pop_hosp_status.shed_col(real_t);
  
  if(output=="simple") {
    
    List L = List::create(_["incidence"] = incidence,
                          _["incidence_intervention"] = incidence_inter,
                          _["susceptible_patients"] = s_patients,
                          _["susceptible_patients_intervention"] = s_patients_inter,
                          _["wards_events"] = ward_event,
                          _["wards_events_intervention"] = ward_event_inter,
                          _["wards_susceptibles"] = ward_event,
                          _["wards_susceptibles_intervention"] = ward_event_inter,
                          _["contaminated_devices"] = cont_mat,
                          _["used_devices"] = used_mat,
                          _["contaminated_devices_intervention"] = cont_mat_inter,
                          _["used_devices_intervention"] = used_mat_inter,
                          _["infection_events_devices"] = inf_eq,
                          _["infection__events_devices_intervention"] = inf_eq_inter,
                          _["count_patients_hospital"] = count_patients,
                          _["newly_admitted_patients"] = new_patient_vec,
                          _["count_admissions_wards"] = count_patient_ward,
                          _["count_tests"]= count_tests,
                          _["device_usage"]= eq_usage );
    
    return L;}
  
}
