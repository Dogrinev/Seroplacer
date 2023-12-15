# Bacterial Serovar Prediction Using the Seroplacer R Package

This repository contains the code for the Seroplacer R package, which makes bacterial serovar predictions using 16S amplicon sequencing data. 

Please see our [vignette](https://dogrinev.github.io/Seroplacer_Vignette.html) for a step-by-step guide on how to use the Seroplacer package with example test data.

## Installation

The easiest way to install Seroplacer is through source code installation directly through Github:

```
library(devtools)
devtools::install_github("dogrinev/Seroplacer")
```

Test data is available pre-built into the package as `Seroplacer::Query_Example` or can be also downloaded as a raw FASTA file from this repository (TEST_DATA.fasta)

## Additional Links

A manuscript describing the development and testing of the Seroplacer package: [Serovar-level Identification of Bacterial Foodborne Pathogens From Full-length 16S rRNA Gene Sequencing](https://pubmed.ncbi.nlm.nih.gov/37425822/)

A [repository](https://github.com/Dogrinev/Seroplacer_Manuscript) containing all reproducible analysis files for the Seroplacer manuscript.
