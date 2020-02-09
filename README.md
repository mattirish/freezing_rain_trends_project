# Scripts & Data for [GLISA](http://glisa.umich.edu/) Great Lakes Freezing Rain Trends Project
This repo contains all the MATLAB and R scripts used in the updated freezing rain/ice storms climatology of the Great Lakes region that I've been working on, along with the data itself.

### Data
- [a_all.mat](a_all.mat): Contains *a*, a struct containing summary data for all 98 observing stations included in the analysis
- [b.mat](b.mat): Contains *b*, a struct containg the full timeseries of freezing rain observations for all 98 stations in the analysis

### Scripts
- [FZRADurationAnalysis.m](FZRADurationAnalysis.m): Analyze spatial and temporal trends in the duration of freezing rain events
- [FZRA_Changepoint_Analysis.R](FZRA_Changepoint_Analysis.R): Perform a changepoint analysis to diagnosis any bias introduced in the transition from manual to automated ASOS freezing rain observations in 1995
- [FZRA_Duration_n_Intensity.m](FZRA_Duration_n_Intensity.m): 
- [HistoricalFZRA_RegionalAvgwProxies.m](HistoricalFZRA_RegionalAvgwProxies.m): bonus! Compare frequency by year with mentions of "freezing rain" in Wikipedia articles on North American ice storms between 2006 and 2010

### Methodology
1. Filter raw ISD observation data to normalize observations to the total number of hours in each year. (IDL scripts by BJ Baule, not included here)
2. Extract all observations for which freezing rain is present. (IDL scripts by BJ Baule, not included here)
3. Do some QC to determine which stations have observations consistent enough for inclusion in the analysis.
	- heatmaps script
	- add spreadsheet with percentage of data available by station
4. Do a changepoint analysis to begin to characterize biases introduced by the 1995 instrumentation change.
	- [FZRA_Changepoint_Analysis.R](FZRA_Changepoint_Analysis.R)
5. Create a struct, *b*, containing all timeseries of freezing rain observations at each station and their metadata, as well as *a*, a struct summarizing station metadata as well as the frequency of freezing rain events and wind conditions.
	- [a_all.mat](a_all.mat)
	- [b.mat](b.mat)
4. Analyze the kinematics of freezing rain—analyze trends in freezing rain frequency and create maps, timeseries plots, etc. 
	- [Matt_spatialtrends_updated.m](Matt_spatialtrends_updated.m)
	- [FZRA_Trend_Analysis_with_modifiedmk.R](FZRA_Trend_Analysis_with_modifiedmk.R)
	- [HistoricalFZRA_RegionalAvgwProxies.m](HistoricalFZRA_RegionalAvgwProxies.m)
5. Create eve

### Software Requirements
- R
	- 
- MATLAB with Mapping Toolbox
	- snctools

*Note:* all third-party R and MATLAB dependencies used in the project can be found at [my other repo here](https://github.com/mattirish/matlab_addons).