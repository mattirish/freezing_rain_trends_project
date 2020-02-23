# Scripts & Data for [GLISA](http://glisa.umich.edu/) Great Lakes Freezing Rain Trends Project
This repo contains all the MATLAB and R scripts used in the updated freezing rain/ice storms climatology of the Great Lakes region that I've been working on, along with the data needed to carry it out.

### Data
- [a_all.mat](a_all.mat): Contains *a*, a struct containing summary data for all 98 observing stations included in the analysis
- [b.mat](b.mat): Contains *b*, a struct containg the full timeseries of freezing rain observations for all 98 stations in the analysis

### Scripts
- [Extraction_All_Stations.m](Extraction_All_Stations.m): Process freezing rain observations from text files in (standard ISD format)[ftp://ftp.ncdc.noaa.gov/pub/data/noaa/ISH-DVD2012/ish-abbreviated.txt] to MATLAB structs *a* and *b* (original ISD text files are too large for inclusion here—contact mairish@umich.edu if you'd like them!)
- [FZRADurationAnalysis.m](FZRADurationAnalysis.m): Analyze spatial and temporal trends in the duration of freezing rain events
- [FZRA_Changepoint_Analysis.R](FZRA_Changepoint_Analysis.R): Perform a changepoint analysis to diagnose any bias introduced in the transition from manual to automated ASOS freezing rain observations in 1995
- [HistoricalFZRA_RegionalAvgwProxies.m](HistoricalFZRA_RegionalAvgwProxies.m): Bonus! Compare frequency by year with mentions of "freezing rain" in Wikipedia articles on North American ice storms between 2006 and 2010
- [FZRA_Intensity_Investigation_CPC.m](FZRA_Intensity_Investigation_CPC.m): 
- [FZRA_NARR_Reanalysis.m](FZRA_NARR_Reanalysis.m): Create maps using NARR (not used for any final figures--this code was adapted elsewhere)

### Methodology
1. Filter raw ISD observation data to normalize observations to the total number of hours in each year and extract all observations for which freezing rain is present. (IDL scripts by BJ Baule, not included here)
2. Process raw text files into structs that can easily be queried.
	- [Extraction_All_Stations.m](Extraction_All_Stations.m)
3. Do some QC to determine which stations have observations consistent enough for inclusion in the analysis.
	- heatmaps script
	- add spreadsheet with percentage of data available by station
4. Do a changepoint analysis to begin to characterize biases introduced by the 1995 instrumentation change.
	- [FZRA_Changepoint_Analysis.R](FZRA_Changepoint_Analysis.R)
5. Create a struct, *b*, containing all timeseries of freezing rain observations at each station and their metadata, as well as *a*, a struct summarizing station metadata as well as the frequency of freezing rain events and wind conditions.
	- [a_all.mat](a_all.mat)
	- [b.mat](b.mat)
6. Analyze trends in the kinematics of freezing rain frequency and create maps, timeseries plots, etc. 
	- [Matt_spatialtrends_updated.m](Matt_spatialtrends_updated.m)
	- [FZRA_Trend_Analysis_with_modifiedmk.R](FZRA_Trend_Analysis_with_modifiedmk.R)
	- [HistoricalFZRA_RegionalAvgwProxies.m](HistoricalFZRA_RegionalAvgwProxies.m)
7. Classify all observations into distinct freezing rain events and record information on their duration and intensity.
	- [FZRA_EventTimes.m](FZRA_EventTimes.m)
	- [FZRA_EventTimes_automated.m](FZRA_EventTimes_automated.m)
8. Analyze trends in the dynamics of freezing rain events and ice storms by pulling synoptic weather data from the NARR reanalysis and carrying out a k-means clustering of synoptic conditions during events.
	- [FZRA_SynopticWeatherTyping_withTwoVariables.m](FZRA_SynopticWeatherTyping_withTwoVariables.m)
	- [FZRA_NARR_Reanalysis.m](FZRA_NARR_Reanalysis.m)

### Software Requirements

- R
	- [R.matlab](https://cran.r-project.org/web/packages/R.matlab/index.html) for transferring data between R and MATLAB
	- [modifiedmk](https://cran.r-project.org/web/packages/modifiedmk/index.html) for Modified Mann-Kendall trend analysis
- MATLAB with Mapping Toolbox
	- [cptcmap](https://github.com/kakearney/cptcmap-pkg/tree/845bf8372ca1be0a83d7c16c05ce3ffaeb033c42)for colormaps
	- [SNCTOOLS](http://mexcdf.sourceforge.net/) for downloading [NARR](https://www.esrl.noaa.gov/psd/data/gridded/data.narr.html) netCDFs via OPeNDAP

All third-party R and MATLAB dependencies used in the project can be found at my other repo [here](https://github.com/mattirish/matlab_addons).