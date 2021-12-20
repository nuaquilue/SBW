# SBW - Spruce budworm outbreaks model

## Introduction

The **SBW** package provides a set of functions to simulate four different phases of spruce budworm outbreaks: pre-epidemic, epidemic, collapse, and calm. It includes landscape-scale processes associated to these phases like neighbour contagion, long-term dispersal, defoliation. The model also incorporates the mortality of forest stands induced by consecutive and severe SBW defoliation seasons.

The SBW model is initialized for the commercial forests of Quebec province (Canada) in 2020. It is a spatially explicit model calibrated to work at 2 km of spatial resolution and 1-year time step.


## Package installation

Users can download and install the latest stable version of the **SBW** package from GitHub as follows (required package devtools should be installed/updated first):

```R
devtools::install_github("nuaquilue/SBW")
```
Additionally, users can have help to run package functions directly as package vignettes, by forcing their inclusion in installation:

```R
devtools::install_github("nuaquilue/SBW", 
                         build_manual = TRUE,
                         build_vignettes = TRUE)
```

## References

Aquilué, N., Bouchard, M., Filotas, É. Endemic, epidemic and collapse phases of spruce budworm outbreaks: a simulation approach (in prep.).
