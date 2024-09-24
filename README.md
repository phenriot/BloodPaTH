
# *BloodPaTH* : Bloodborne Pathogens Transmission in Hospitals 

## Introduction <a href="README.md"> <img src="https://github.com/phenriot/BloodPaTH/blob/main/Other/BloodPaTH_logo.png" align="right" width="150"/> </a>

The *BloodPaTh* package is a C++ coded tool allowing the investigation of bloodborne pathogens transmission within healthcare settings using longitudinal prospective data. This agent-based model reproduces : </br>
* **i.** geographical movements of patients between wards within the hospital, </br> 
* **ii.** the medical devices dynamics (contamination, use and reuse of devices as well as the variation of sterile devices availability), and </br>
* **iii.** the epidemic dynamics of transmission for patients undergoing different types of invasive procedures (following a SEI model). </br> </br>

More details about the model are available in [*An agent-based model to simulate the transmission dynamics of bloodborne pathogens within hospitals*](https://www.medrxiv.org/content/10.1101/2023.11.14.23298506v1).

## Installation 

The *BloodPaTH* package can be installed in two ways : 

* **1.** Simply run the following command line within your R session :
``` r
devtools::install_github("phenriot/BloodPaTH")
```

* **2.** After downloading the BloodPaTH_1.0.tar.gz file, install it in your R environment using the following command line:
 ``` r
 tools::Rcmd("INSTALL BloodPaTH_1.0.tar.gz")
 ```

**NB** : You will need te following packages to proceed : *tools*, *truncnorm*, *mc2d*, *RcppArmadillo*, *Rcpp*

## How does it work ? 

### List of model parameters
* ***time_step*** : time-step (must be expressed in hours; **example** : if time step is 1 minute, inform 1/60) | <ins> type = float </ins>
* ***t*** : simulation time (must be expressed in hours; **example** : if simulation time is 1 year, inform 24*365 = 8,760) | <ins> type = float </ins>
* ***pathogen*** : bloodborne pathogen type (must be "HCV", "HBV" or "custom"; if "custom" then the parameter *dist_risk* needs to be informed; The "HIV" option is coming soon) | <ins> type = string </ins>
* ***min_e_phase*** : minimum value of the [eclipse phase](https://bio.libretexts.org/Courses/Portland_Community_College/Cascade_Microbiology/06%3A_Acellular_Pathogens/6.2%3A_The_Viral_Life_Cycle) (must be expressed in hours) | <ins> type = float </ins>
* ***max_e_phase*** : maximum value of the [eclipse phase](https://bio.libretexts.org/Courses/Portland_Community_College/Cascade_Microbiology/06%3A_Acellular_Pathogens/6.2%3A_The_Viral_Life_Cycle) (must be expressed in hours) | <ins> type = float </ins> 
* ***nb_patients*** : number of patients | <ins> type = integer </ins>
* ***nb_wards*** : number of wards within the healthcare setting | <ins> type = integer </ins>
* ***nb_adm*** : number of admission routes (i.e, most of the time the number of departments) | <ins> type = integer </ins>
* ***adm_prob*** : probability of admission in each of the admission route (must be the same size as ***nb_adm*** ; this vector needs to sum to 1) | <ins> type = float vector </ins>
* ***prev_type*** : you can either inform a prevalence at the admission route level ("admission route") or at a ward level ("ward") | <ins> type = string </ins>
* ***prev_init*** : prevalences upon admission (if ***prev_type*** = "admission route", it needs to be the same size as ***nb_adm*** ; if ***prev_type*** = "ward", it needs to be the same size as ***nb_wards***) | <ins> type = float vector </ins>
* ***WT_matrix*** : merged transition matrices for all admission routes; the *list_to_combined_matrices* R function helps to convert a list of transition matrices as follows &rarr; `list_to_combined_matrices(input = your_list, type = 1, nb_wards = nb_wards)` . Transition matrices must be of size (***nb_wards***+1) * (***nb_wards***+1) (i.e, include an extra row and an extra column for the "discharged" event)  | <ins> type = float matrix </ins>
* ***nb_procedures*** : number of procedure types performed within the healthcare setting (must include the "no procedure event"; **example** : if there are 10 different types of procedures perfomed within the hospital then inform 10+1 = 11) | <ins> type = integer </ins>
* ***procedure_names*** : vector of procedure names | <ins> type = string vector </ins>
* ***nb_devices*** : number of different device types used within the healthcare setting | <ins> type = integer </ins>
* ***device_names*** : vector of devices names | <ins> type = string vector </ins>
* ***nb_devices_new*** : initial number of new sterile devices in each ward; matrix of size ***nb_wards*** * ***nb_devices***  | <ins> type = integer matrix </ins>
* ***nb_devices_used*** : initial number of non-sterile devices (previously used) in each ward; matrix of size ***nb_wards*** * ***nb_devices***  | <ins> type = integer matrix </ins>
* ***nb_devices_cont*** : initial number of non-sterile and contaminated devices in each ward; matrix of size ***nb_wards*** * ***nb_devices***  | <ins> type = integer matrix </ins>
* ***refill_quantities*** : quantity of devices to add to the pool of sterile devices for each type at each refill event; matrix of size ***nb_wards*** * ***nb_devices***  | <ins> type = integer matrix </ins>
* ***refill_freq*** : refill frequency for each type of device in each ward (must be expressed in hours); matrix of size ***nb_wards*** * ***nb_devices***  | <ins> type = integer matrix </ins>
* ***sterilization_prob*** : efficient sterilization probability for each type of device (might evolve to a matrix of ward-specific probabilities in the future) | <ins> type = float vector </ins> 
* ***PPM_matrix*** : merged matrices of probabilities of undergoing different procedures (columns) within each ward (rows);  the *list_to_combined_matrices* R function helps to convert a list of probability matrices as follows &rarr; `list_to_combined_matrices(input = your_list, type = 2, nb_procedures = nb_procedures)` . These probability matrices must be of size ***nb_wards*** * ***nb_procedures*** | <ins> type = float matrix </ins>
* ***dist_risk*** : parameters of the distributions of the of risk of getting infected for each of the procedures. The first three columns correspond to the values of the parameters of these distributions and the fourth to the name of these distributions (log-normal : "lnorm", gaussian :"norm" or PERT : "pert"). Only the two first columns are used for the [log-normal](https://en.wikipedia.org/wiki/Log-normal_distribution) and [normal](https://en.wikipedia.org/wiki/Normal_distribution) distributions (for $\mu$ and  	$\sigma^2$); all columns are used for the [PERT](https://en.wikipedia.org/wiki/PERT_distribution) distribution (for $a$, $b$ and $c$). This data frame needs to have ***nb_procedures*** rows; the last corresponds to the risk of getting infected when undergoing no procedure (which is 0), then it can be defined as something like &rarr; `rbind(your_dist_risk_data_frame,c(0,0,NA,"lnorm"))` | <ins> type = data frame </ins>
* ***intervention*** : should an intervention occur to reduce the burden of infections ? Can either be "none", "upon-admission" or "ward-level". The upon-admission intervention corresponds to a systematic screening of patients upon-admission at the healthcare setting level with a probability of ***prob_screening***. The ward-level intervention corresponds to a systematic screening of patients entering for the first time the most at risk wards (defined with the parameter ***at_risk_wards***) with a probabilty of ***prob_screening***. For each case, the probability of efficient sterilization for devices used on positive screened patients is set to 1.  | <ins> type = string </ins>
* ***prob_screening*** : probability of patient screening upon-admission at a healthcare setting or ward level. | <ins> type = float </ins>
* ***at_risk_wards*** : vector of the most at risk wards (i.e, wards in which we want to set an intervention); the wards needs to be informed as integers, corresponding to their rank (or ID). | <ins> type = integer vector </ins>
* ***output*** : should the output be "simple" or "detailed" (the second option is not implemented yet). | <ins> type = string </ins>
* ***id_sim*** : simulation ID; this is used to facilitate replication of results and to make baseline and intervention scenarios comparable. If you are running multiple simulations, you will need to inform this parameter with the simulation number to allow replication.
### Model output (simple) 

*... Coming soon ...*

### Model output (detailed)

*... Coming soon ...*

**NB**: This package is an alpha version and is still under development. A Shiny app is currently being developped. 

## Application

You will find an application of the model within the *Application* folder. 

The *model_application_example.R* file allows you to understand what are the parameters of the model and to use the model to investigate the spread of HCV in a synthetic healthcare setting (*see below*).  

To run the code you will need to download the *Data example* folder, in which you will find synthetic data :
- *List_Transition_Matrices.rds* : An RDS file containing a list of 2 transition matrices between wards of size 29x29
- *List_Proc_Prob_Matrices.rds* : An RDS file containing a list of 2 matrices of probabilites of undergoing a set of procedures while being hospitalized in each of the wards (size : 28x10)
- *association_devices_procedures.csv* : A CSV file summarising the association between devices (column 'ID_devices') and procedures (column 'ID_procedures')
- *risk_dist.csv* : A CSV file detailing the parameters of the distribution of the risk of getting infected for each type of procedure 

**NB:** You will have to change the path when loading the data within you R session.  

*Analyses_PlosCompBiol_paper.R* is the R file summarising all the analyses performed to obtain the results presented in [*An agent-based model to simulate the transmission dynamics of bloodborne pathogens within hospitals*](https://www.medrxiv.org/content/10.1101/2023.11.14.23298506v1). Unfortunately, it cannot be used in its current form because the data used in this project is considered sensitive and cannot be shared publicly. 

## Contact

If you have any questions, please reach the author <a href="mailto:paul.henriot@protonmail.com">Paul Henriot</a>
