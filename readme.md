# Identifying `useful' fitness models: balancing the benefits of added complexity with realistic data requirements in models of individual plant fitness
## Code and data used in publication.
#### Trace Martyn*, Daniel Stouffer, Oscar Godoy, Ignasi Bartomeus, Abigail Pastore, Margaret Mayfield
#### *Corresponding author; E-mail: martyn.ecology@gmail.com
##### Initial commit: April 22, 2019
##### Most recent commit: Dec 9, 2020

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
