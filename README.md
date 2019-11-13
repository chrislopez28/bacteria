(Documentation in progress)

# Bacteria - Functions for Indicator Bacteria Compliance Analysis in California

Bacteria is an R package that contains functions that can be used to evaluate compliance with water quality objectives for indicator bacteria (_E. coli_, _Enterococcus_, Total Coliform, and Fecal Coliform). The primary use for these functions is for the evaluation of bacteria water quality monitoring data within the jursidction of the Los Angeles Regional Water Quality Control Board. 

## Installation

1. If not installed already, install the `devtools` package. This can be done using the following: `install.packages("devtools")`
2. Install the `bacteria` package with the following: `devtools::install_github("chrislopez28/bacteria")`. 

## Using the `bacteria` package

### bact_check

The `bact_check` function takes in a data frame of bacteria monitoring data and outputs a list of data frames (each element of the list corresponding to a unique monitoring location) that contains the results of the analysis.

#### Format of input `df`
The `bact_check` requires that the initial input data frame has the following variables:

* StationCode - a character vector
* SampleDate - a date vector
* WeatherCondition - a character vector
* AnalyteName - a character vector
* Result - a numeric vector
* ResQualCode - a character vector
* MDL - a numeric vector
* RL - a numeric vector

### ann_exceeds

The `ann_exceeds` function calculates the number of _Annual Exceedance Days_. The arguments for this function are the same as that of the `bact_check` function.

### Temporal Heatmap Visualizations

The `bact_heatmap` function creates a heatmap that summarizes compliance with bacteria water quality objectives while also providing information on sampling dates and weather conditions. 

![An Example Temporal Heatmap](https://github.com/chrislopez28/bacteria/blob/master/heatmap5.png "Example Heatmap")
