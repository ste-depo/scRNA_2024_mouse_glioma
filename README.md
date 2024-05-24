# scRNA analysis of treated mouse glioma samples

Code, data and figures that reproduce the analyses present in the study:

_XXX_

* Pipeline.Rmd:
  Script for importing data from cellranger output, perform filtering and integration

* Pipeline.html:
  Knit to HTML output of Pipeline.Rmd 

* SourceFunctions.R
  Custom R functions that are required to execute Pipeline.Rmd

* run_souporcell.sh:
  Code for souporcell analysis

* InputData:
  cellranger outputs (empty, they should be filled GSE256358 processed data) and souporcell outputs

* Figures:
  figures that the script can produce
