# tamu_daily_lake_evap
The code provided here is for daily lake evaporation estimations, and is based on the algorithm described in the following paper: 
 * Zhao, G., & Gao, H. (2019). Estimating reservoir evaporation losses for the United States: Fusing remote sensing and 
   modeling approaches. Remote sensing of Environment, 226, 109-124. 
 * Zhao, B., Huntington, J. Pearson, C., Zhao, G., Ott, T., ..., Gao, H. Developing a General Daily 	Lake Evaporation Model 
   and Demonstrating its Application in the State of Texas. Water Resource Research. (2024).

Usage Restriction:
* The code provided here is intended for academic and research purposes only.  When using this code, please cite the aforementioned papers. It should not be redistributed or used for commercial purposes without explicit permission from the original authors. 

Tips:
* Please execute the main function in 'DLEM_run.py' to operate the daily lake evaporation model. 
* We have included an example using Lake Powell, which requires three essential inputs: 'powell_area_depth.csv' describing the dynamic area and depth; 'powell_meteo.csv' detailing the primary meteorological forcings extracted from GridMET and RTMA; and 'powell_fetch.csv' providing the fetch length in eight different directions. 
* The file 'powell_DLEM_output.csv' will contain the DLEM evaporation estimates.
