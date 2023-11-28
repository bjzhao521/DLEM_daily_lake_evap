# tamu_daily_lake_evap
The code for daily lake evaporation estimations.
 
This code is based on the algorithm described in the paper: 
 * Zhao, G., & Gao, H. (2019). Estimating reservoir evaporation losses for the United States: Fusing remote sensing and 
 * modeling approaches. Remote sensing of Environment, 226, 109-124. 
 * Zhao, B., Huntington, J. Pearson, C., Zhao, G., Ott, T., ..., Gao, H. Developing a General Daily 	Lake Evaporation Model 
 * and Demonstrating its Application in the State of Texas. Water 	Resource Research. (2023b), In review.

Usage Restriction:
 * 
 * The code provided here is intended for academic and research purposes only. You should cite the aforementioned papers when 
 * using this code. It is not to be redistributed or used for commercial purposes without explicit permission from the 
 * original authors. 

Tips:
* 1. Please execute the main function in 'DLEM_run.py' to operate the daily lake evaporation model. 
* 2. We have included an example using Lake Powell, which requires three essential inputs: 'powell_area_depth.csv' describing * the dynamic area and depth; 'powell_meteo.csv' detailing the primary meteorological forcings extracted from GridMET and     * RTMA; and 'powell_fetch.csv' providing the fetch length in eight different directions. 
* 3. The file 'powell_DLEM_output.csv' will contain the DLEM evaporation estimates.
