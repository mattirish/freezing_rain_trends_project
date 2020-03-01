# Changepoint Analysis for ISD Historical Freezing Rain Time Series
# Matt Irish
# July 2019

library(ecp)
library(R.matlab)
library(ggplot2)

setwd('/Users/mattirish/Documents/MATLAB/freezing_rain_trends_project')
a_all <- readMat("a_all.mat")
YearFreq_rel <- a_all[["a"]][[11]]
YearFreq_rel_domain_avg <- apply(YearFreq_rel,2, median, na.rm = T)
MonthFreq_yearly_rel_series <-  a_all[["a"]][[12]]
MonthFreq_yearly_rel_series_domain_avg <- apply(MonthFreq_yearly_rel_series,2, median, na.rm = T) #Domain average freezing rain for every month in series
MonthFreq_yearly_rel_series_domain_sum <- apply(MonthFreq_yearly_rel_series,2, sum, na.rm = T) #Domain average freezing rain for every month in series

# Do a changepoint detection on the yearly frequencies across the region and at each station:
avg_changepoints <- e.divisive(as.matrix(colMeans(YearFreq_rel,na.rm= T)), sig.lvl=.999, R=199, k=NULL, min.size=11, alpha=1)

## Plot the domain average and means across the changepoint:
end_yr_of_first_period <- seq(1976,2014)[avg_changepoints$estimates[2]-1]
period1_avg <- mean(colMeans(YearFreq_rel,na.rm= T)[1:avg_changepoints$estimates[2]-1])
period2_avg <- mean(colMeans(YearFreq_rel,na.rm= T)[avg_changepoints$estimates[2]:39])
sprintf('End year of first period is %s and p = %s.',end_yr_of_first_period,avg_changepoints$p.values)
sprintf('Period 1 avg: %s',period1_avg)
sprintf('Period 2 avg: %s',period2_avg)
plot(data.frame(seq(1976,2014),colMeans(YearFreq_rel,na.rm= T)),type='l',xlab="Year",ylab="Hours of Freezing Rain per Year")
lines(data.frame(seq(1976,end_yr_of_first_period),rep(period1_avg,avg_changepoints$estimates[2]-1)),col="blue")
lines(data.frame(seq(end_yr_of_first_period+1,2014),rep(period2_avg,40-avg_changepoints$estimates[2])),col="red")

avg_changepoints_monthly <- e.divisive(as.matrix(colMeans(MonthFreq_yearly_rel_series,na.rm= T)), sig.lvl=.8, R=199, k=NULL, min.size=48, alpha=1)

#changepoint_onestationtest <- e.divisive(as.matrix(na.omit(YearFreq_rel[1,])), sig.lvl=.04, R=199, k=NULL, min.size=2, alpha=1)

changepoints <- apply(YearFreq_rel, 1, function(x) e.divisive(as.matrix(na.omit(x)), sig.lvl=.05, R=199, k=NULL, min.size=10, alpha=1))

changepoints_dt <- as.data.table(do.call(rbind, changepoints))

changepoint_year_indices <- unlist(changepoints_dt[,estimates])
changepoint_year_indices <- changepoint_year_indices[! changepoint_year_indices %in% c(1,40)]
p <- hist(changepoint_year_indices,breaks=diff(range(changepoint_year_indices)))


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