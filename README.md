
# *BloodPaTH* : Bloodborne Pathogens Transmission in Hospitals 

## Introduction <a href="README.md"> <img src="https://github.com/phenriot/BloodPaTH/blob/main/Other/BloodPaTH_logo.png" align="right" width="150"/> </a>

The *BloodPaTh* package is a tool allowing the investigation of bloodborne pathogens transmission within healthcare settings using longitudinal prospective data.

## Installation 

After downloading the *BloodPaTH* package (BloodPaTH_1.0.tar.gz) please start by installing it in your R environment using the following command line (you will need the *tools* package to process):
 > tools::Rcmd("INSTALL BloodPaTH_1.0.tar.gz")

## How does it work ? 

The main model was coded in C++.

*... Coming soon ...*

**NB**: This package is in alpha version and is still under development. A Shiny app is currently being developped. 

## Application

You will find an application of the model in the *model_application_example.R* file. 

To run the code you will need to download the *Data example* folder, in which you will find synthetic data :
- *List_Transition_Matrices.rds* : An RDS file containng a list of 2 transition matrices between wards of size 29x29
- *List_Proc_Prob_Matrices.rds* : An RDS file containng a list of 2 matrices of probabilites of undergoing a set of procedures while being hospitalized in each of the wards (size : 28x10)
- *association_devices_procedures.csv* : A CSV file summarising the association between devices (column 'ID_devices') and procedures (column 'ID_procedures')
- *risk_dist.csv* : A CSV file detailing the parameters of the distribution of the risk of getting infected for each type of procedure 

**NB:** You will have to change the path when loading the data within you R session. 

You will find in the *Functions_plot.R* file the functions allowing visualization of the output. 

## Other info 

*Analyses_PlosCompBiol_paper.R* is the R file summarising all the analyses performed to obtain the results presented in [*An agent-based model to simulate the transmission dynamics of bloodborne pathogens within hospitals*](https://www.medrxiv.org/content/10.1101/2023.11.14.23298506v1).
Unfortunately, it cannot be used in its current form because the data used in this project is considered sensitive and cannot be shared publicly. 

If you have any question, please reach the author <a href="mailto:paul.henriot@protonmail.com">Paul Henriot</a>
