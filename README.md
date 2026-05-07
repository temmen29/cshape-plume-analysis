# cshape-plume-analysis
Code and analysis for entraining plume simulations using tropical ARMBE and COSMIC-2 data
# cshape-plume-analysis

Code and analysis for entraining plume simulations using tropical ARM and COSMIC-2 data.

This repository contains processing scripts, plume-model calculations, and analysis notebooks associated with:

> *The C-Shape in Tropical Thermodynamic Profiles: An Explanation Combining Entrainment and Condensate Loss*

---

## Repository contents

### Core functions and model

- `beam_model.py`  
  Bulk Entropy Accounting Model (BEAM) used for entraining plume calculations and condensate-loss experiments. Includes moist thermodynamic utilities, equivalent potential temperature calculations, mixed-phase partitioning, condensate retention parameterizations. 

---

### Processing scripts

- `ARM_process_theta_e.py`  
  Processes ARM/ARMBE at Manus and Nauru site thermodynamic profiles and condensate from ARM MICROBASE. Computes variables used in the plume analysis, including equivalent potential temperature, saturation equivalent potential temperature, humidity, condensate, etc. Writes netCDF for workflow.

- `cosmic_2_process_script.py`  
  Processes COSMIC-2 radio occultation profiles and interpolates thermodynamic variables onto a common pressure grid. Writes netCDF for workflow.

- `imerg_match_to_COSMIC2.py`  
  Matches COSMIC-2 profiles to nearby IMERG precipitation. Writes netCDF for workflow.

---

## Analysis notebooks

- `condensate_loss_mixing_ARMBE.ipynb`  
  ARM entraining plume analysis, including condensate loss- entrainment compensation experiments with lower-free-tropospheric buoyancy constraints.

- `condensate_loss_mixing_COSMIC2.ipynb`  
  Analysis of COSMIC-2 thermodynamic structure and precipitation-conditioned $\theta_{es}$ anomalies.

---

## Data sources

Raw datasets are publicly available from their original archives:

- ARM Best Estimate (ARMBE): https://www.arm.gov/data/science-data-products/vaps/armbe
- ARM MICROBASE: https://adc.arm.gov/discovery/#/results/s::microbase
- COSMIC-2 radio occultation: https://doi.org/10.5065/t353-c093
- GPM IMERG precipitation: https://gpm.nasa.gov/data

Raw datasets are not redistributed here.

Processed intermediate NetCDF files may be generated using the included scripts.

---

## Suggested workflow

### 1. 'ARM_process_theta_e.py' -- Process ARM thermodynamic and condensate variables into netCDF file
### 2. 'cosmic_2_process_script.py' -- Process COSMIC-2 radio occultation data into netCDF file
### 3. 'imerg_match_to_COSMIC2.py' -- Match IMERG precipitation measurements to COSMIC-2 ROs
### 4. 'condensate_loss_mixing_ARMBE.ipynb' -- Examples of running the BEAM model using ARMBE data and analysis
### 5. 'condensate_loss_mixing_COSMIC2.ipynb' -- Example of analysis of $\theta_e$ based analysis using COSMIC-2 matched with IMERG 
