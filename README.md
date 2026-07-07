# Integration-of-metabolomics-into-early-phase-clinical-trials-from-2010-to-2025

This repository contains the R code used for the analysis reported in:

Integration of metabolomics into early phase clinical trials from 2010 to 2025

The analysis identifies early phase interventional clinical trials registered in the WHO International Clinical Trials Registry Platform (ICTRP), classifies trials by phase, registry, country, year and outcome field, and identifies studies incorporating metabolomics or related omics terminology.

Repository contents

•	analysis.R
Main R script used to clean the ICTRP dataset, identify early phase trials, detect metabolomics-related studies, remove duplicate registrations, generate summary tables and produce figures.

•	country_lookup.csv
Lookup table used to standardise inconsistent country names in ICTRP records.

•	continent_analysis.csv
Manually curated country-to-continent classification and continent-level summary input.

•	outputs/
Folder for generated analysis outputs, including cleaned datasets, summary tables and figures.

Data source

The analysis uses the WHO ICTRP clinical trial dataset exported as:

ICTRP_Full_R.csv

This file is not included in the repository because it should be obtained directly from the WHO ICTRP or relevant permitted source. To reproduce the analysis, place ICTRP_Full_R.csv in the project root directory before running the script.

The input dataset used in the original analysis contained 1,173,018 records and 63 variables.

Main outputs

The script generates:

•	ICTRP_All.csv
Early phase interventional clinical trials from 2010 to 2025 after duplicate handling.

•	ICTRP_Met.csv
Subset of early phase trials containing metabolomics-related search terms.

•	country_met_fraction.csv
Country-level counts and proportions of metabolomics studies.


Additional outputs include summary tables and figures stratified by:
•	country
•	continent
•	clinical trial phase
•	registry
•	year
•	outcome classification

Search strategy

Metabolomics-related studies were identified by searching selected ICTRP fields for the following terms:
•	metabolom*
•	metabonom*
•	metabolite profil*
•	metabolite fingerprint*
•	small molecule profiling
•	pharmacometabolom*
•	lipidom*
•	glycom*
•	steroidom*
•	metaboli* phenotyping

The searched fields were:
•	public title
•	scientific title
•	primary outcome
•	secondary outcomes
•	conditions
•	interventions

Software requirements
The analysis was conducted in R using the following packages:
•	tidyverse
•	dplyr
•	ggplot2
•	stringr
•	readr
•	janitor
•	lubridate
•	purrr
•	tidyr
•	tibble
•	ggrepel
Package versions should be recorded using sessionInfo() or an renv.lock file for full reproducibility.

How to run the analysis
1.	Clone or download this repository.
2.	Place ICTRP_Full_R.csv in the project root directory.
3.	Ensure country_lookup.csv and continent_analysis.csv are present.
4.	Open the R project or set the working directory to the repository root.
5.	Run analysis.R.
Example:
source("analysis.R")

Reproducibility notes
The ICTRP country field contains heterogeneous free-text entries. Country names were therefore standardised using a manually curated lookup table. Trials hosted in multiple countries were handled using fractional weighting so that each trial contributed a total weight of one across hosting countries.
Duplicate registrations were identified using normalised public and scientific titles, with the earliest registration date retained where duplicate records were detected.

Citation
If using or adapting this code, please cite the associated manuscript and the archived repository DOI.
Thompson J, et al. Integration of metabolomics into early phase clinical trials from 2010 to 2025.

Licence
This code is shared under the MIT License, unless otherwise stated.

Contact
For questions about the analysis, please contact the corresponding author of the associated manuscript.


