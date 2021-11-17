# Bayesian Data Analysis Project

_Yann McLatchie, Arina Odnoblyudova_

A project completed for the Bayesian Data Analysis module at Aalto University (2021).

## How to compile the project

The project is organised into some subdirectories:

- `R`: Contains R scripts for general data analysis and model building
- `stan`: Contains Stan files for use in the project
- `rmd`: Contrains the Rmd files for sections of the final project

In order to knit the complete report, write Rmd code into the child files found in `rmd`, embed them into `lung_survival.Rmd`, and then knit the parent `lung_survival.Rmd` which produces the final PDF output.
