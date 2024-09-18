
# *BloodPaTH* : Bloodborne Pathogens Transmission in Hospitals 

## Introduction <a href="README.md"> <img src="https://github.com/phenriot/BloodPaTH/blob/main/Other/BloodPaTH_logo.png" align="right" width="150"/> </a>

The *BloodPaTh* package is a tool allowing the investigation of bloodborne pathogens transmission within healthcare settings using longitudinal prospective data. This agent-based model reproduces : </br>
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

The main model was coded in C++.

*... More coming soon ...*

**NB**: This package is an alpha version and is still under development. A Shiny app is currently being developped. 

## Application

You will find an application of the model within the *Application* folder. 

The *model_application_example.R* file allows you to understand what are the parameters of the model and to use the model to investigate HCV spread in a synthetic healthcare setting (*see below*).  

To run the code you will need to download the *Data example* folder, in which you will find synthetic data :
- *List_Transition_Matrices.rds* : An RDS file containng a list of 2 transition matrices between wards of size 29x29
- *List_Proc_Prob_Matrices.rds* : An RDS file containng a list of 2 matrices of probabilites of undergoing a set of procedures while being hospitalized in each of the wards (size : 28x10)
- *association_devices_procedures.csv* : A CSV file summarising the association between devices (column 'ID_devices') and procedures (column 'ID_procedures')
- *risk_dist.csv* : A CSV file detailing the parameters of the distribution of the risk of getting infected for each type of procedure 

**NB:** You will have to change the path when loading the data within you R session.  

*Analyses_PlosCompBiol_paper.R* is the R file summarising all the analyses performed to obtain the results presented in [*An agent-based model to simulate the transmission dynamics of bloodborne pathogens within hospitals*](https://www.medrxiv.org/content/10.1101/2023.11.14.23298506v1). Unfortunately, it cannot be used in its current form because the data used in this project is considered sensitive and cannot be shared publicly. 

## Contact

If you have any question, please reach the author <a href="mailto:paul.henriot@protonmail.com">Paul Henriot</a>
