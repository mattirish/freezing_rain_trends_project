# Trend Analysis for ISD Historical Freezing Rain Time Series
# Matt Irish
# May 2019

library(modifiedmk)
library(R.matlab)
library(ggplot2)

setwd('/Users/mattirish/Documents/MATLAB/freezing_rain_trends_project')
a_all <- readMat("a_all.mat")
YearFreq_rel <- a_all[["a"]][[11]]
YearFreq_rel_domain_avg <- apply(YearFreq_rel,2, median, na.rm = T)
MonthFreq_yearly_rel_series <-  a_all[["a"]][[12]]
MonthFreq_yearly_rel_series_domain_avg <- apply(MonthFreq_yearly_rel_series,2, median, na.rm = T) #Domain average freezing rain for every month in series
MonthFreq_yearly_rel_series_domain_sum <- apply(MonthFreq_yearly_rel_series,2, sum, na.rm = T) #Domain average freezing rain for every month in series


# Estimate the annual trend for each station:
yearly_trends <-apply(YearFreq_rel, 1, function(x) mmkh(x,c=0.95))
yearly_trends_3lag <-apply(YearFreq_rel, 1, function(x) mmkh3lag(x,c=0.95))
# Estimate the annual trend for the entire region:
yearly_trends_domain_avg <- mmkh(YearFreq_rel_domain_avg,c=0.95)

writeMat("yearly_trends.mat",yearly_trends_data=yearly_trends, yearly_trends_3lag_data=yearly_trends_3lag,fixNames=T)


# Estimate the monthly trend for each station:
pval_monthly_domain_avg <- vector("numeric", 12)
sen_monthly_domain_avg <- vector("numeric", 12)
for (month in 1:5){
  monthly_trends_domain_sum <- mmkh(MonthFreq_yearly_rel_series_domain_avg[seq(month,ncol(MonthFreq_yearly_rel_series), by=12)], c=0.95)
  
  pval_monthly_domain_avg[month] <- monthly_trends_domain_sum[2]
  sen_monthly_domain_avg[month] <- monthly_trends_domain_sum[7]
}


ggplot(yearly_trends[2,], aes(1:length(yearly_trends[2,]), yearly_trends[2,] )) +
  ggline()


## DIY-MODIFIED SEASONAL (HAMED et. al):
# Run on each month separately
# Add Z stats together and var(Z)s
# Compute new overall Z stat for each station
# Compare with standard normal distribution for significance.

Zc <- vector("numeric", 12)
sen <- vector("numeric", 12)
varZc <- 
yearly_trends_seasonal <- vector("numeric", 12)
yearly_trends_seasonal <- vector("numeric", 12)
for (m in 1:97){
for (month in 1:12){
  yearly_trends_seasonal[month] <- mmkh(MonthFreq_yearly_rel_series[m,seq(month,ncol(MonthFreq_yearly_rel_series), by=12)], c=0.95)

}
}