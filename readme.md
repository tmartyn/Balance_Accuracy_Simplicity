# Identifying `useful' fitness models: balancing the benefits of added complexity with realistic data requirements in models of individual plant fitness
## Code and data used in publication.
#### Trace Martyn*, Daniel Stouffer, Oscar Godoy, Ignasi Bartomeus, Abigail Pastore, Margaret Mayfield
#### *Corresponding author; E-mail: martyn.ecology@gmail.com
##### Initial commit: April 22, 2019
##### Most recent commit: Dec 18, 2020

## Introduction

This is a project trying to understand how best to balance model complexityand data tratibility in models that include higher-order interactions.

## Repository structure

Present are five folders:
 - Focal.ID - folder that holds code to run models for focal-species identity analysis
 - Neigh.ID - folder that holds code to run models for neighbor-species identity analysis
 - Species.grouping - folder that holds code to group neighbor species by rarity and functional traits
 - figures - folder that holds output figures
 - raw.data - folder that holds raw data and cleaning code

## To run code

To run the analyses, please open the "FINAL.Rproj" file which will open R.Studio and set the working directory from the folder that contains the project file. Then open the "run.analyses.R" file. From this file all data cleaning and analysis code is sourced and run - running this file should run all analyses. If wanting to run the analyses for seperate focal species, use the code with the 'by.focal' naming.

## Data 

Data is accessible at DOI: https://doi.org/10.5061/dryad.zs7h44j7f.

### Data headers

#### mayfield.Rdata:
 - List of dataframe for each of the six focal species
	- Rows are individual plants and columns are:
		- "focal" (the name of the focal species)
		- "seeds" (the number of seeds that individual produced)
		- "quadrat" (the quadrat label) 
		- "site" (the site label - either K for Kunjin or B for Bendering)
		- a list of columns that are headed with the names of neighboring species and include the abundance of that species.

#### Caracoles.competition.update.14.05.2018.csv
 - Rows are individual plants and columns include:
  - "order" (the order of rows that the data was input into the spreadsheet)
  - "date" (the date the data was input into the spreadsheet)
  - "plot" (the plot number 1-8)
  - "subplot" (the subplot code A-F, 1-6)
  - "focal" (the focal species code, see trait data for species name)
  - "fruit" (the number of fruit for the individual plant)
  - "seed"	(the number of seed for the individual plant)
  - a list of columns that are headed with the names of neighboring species and include the abundance of that species.

#### species.traits.Australia.Spain.2019.csv
The two datasets, Spain and Australia were analyzed seperately from one another, therefore we decided to keep the data as collected by the respective groups - therefore units may differ between the two but noted below)
 - Rows are different species and columns include:
	- "Data.name" (name given to the the data based on researcher)
	- "Data" (name given to the data based on location)
	- "SP.code" (four letter code given to the species)
	- "Species" (species scientific name)
	- "native_exotic" (whether the species is native or exotic to the area)
	- "grass_forb" (whether the species is a grass or forb)
	- "annual_perennial" (whether the species is an annual or perennial)
	- "family" (the plant family the species belongs to)
	- "species.dataset.code" (the code used for that species in the neighborhood datasets)
	- "exotic" (whether the species is native or exotic to the area, 1 = exotic, 0 = native)
	- "measured.sla" (mean sla in cm^2/g for Spain and mm^2/mg for Australia)
	- "mean.seed.mass" (mean seed mass in grams)
	- "measured.height" (mean height of species in cm for Spain and mm for Australia)
	- "measured.width" (mean longest radius in mm)
	- "leaf.area.for.sla" (mean leaf area cm^2 for Spain and mm^2 for Australia)
	- "leaf.dry.mass" (mean leaf dry mass in g)
	- "leaves.trait.handbook" (qualitative leaf trait based off Perez-Harguindeguy et al (2013) New handbook for standardize measurement of plant functional traits worldwide. Australian Journal of Botnay 61, 167-234)
	- "vertical_structure" (ratio of height to width of the speices, closer to 1.0 indicates taller than wide)
	- "leaf_area_index" (ratio of leaf area to ground area coverage)
	- "root_diamter_at_base_of_the_stem" (root diameter in cm)
	- "specific_root_length" (ratio of root length to dry mass of fine roots)
	- "volumetric_root_density" (total length of roots per unit soil volume)
	- "specific_root_area" (ratio of root area to dry mass of fine roots)
	- "carbon/nitrogen" (carbon nitrogen ratio of plant)
	- "C13" (the fractionation of isotope C_13)
	- "N15" (the fractionation of isotope N_15)


[![DOI](https://zenodo.org/badge/182592357.svg)](https://zenodo.org/badge/latestdoi/182592357)
