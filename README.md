"# Gola_CocoaBugs" 

Scripts are ordered 1-9

1: Imports and extracts required insect data from raw excel files provide

2: Imports coordinates of trap points and various satellite rasters
  - Calculates mean/sd pixel values for a range of buffer distances around each trap
  - Creates GLCM rasters & calculates statistics for each trap
  - Creates canonical components from above values
  
3: Fits and cross validates a joint species distribution model using environmental data measured in the field

4: Fits a joint species distribution model using satellite derived variables

5: Cross validates jsdm fitted with satellite data

6: Create satellite variables required to predict across the GRNP and leakage belt using the satellite data model

7: Predict species level community composition across the GRNP and leakage belt, then calculate & plot predicted species richness

8: Using species predictions, calculate and map similarity/differences in community composition across the GRNP and leakage belt

9: Predict GEDI measured elements of forest structure using the same variables the joint species distribution model was fit with
  - Note: These are different from the trap site points & predictor values were calculated seperately




DATA
The data2 folder includes all input and output data used and created in running the code
  - NB: it excludes satellite raster inputs due to size constraints, but includes the calculated outputs for the next stages


NON-CRAN PACKAGES
Requires sgdm to create canonical components
  - https://github.com/sparsegdm/sgdm_package

Requires sjSDM package, which can be loaded from CRAN, but the github repository provides important information
  - https://github.com/TheoreticalEcology/s-jSDM

