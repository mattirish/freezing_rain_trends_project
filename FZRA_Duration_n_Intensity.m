%FZRA_Duration_n_Intensity

load a_all
load b
stationnames = fieldnames(b);

%Import station regional cluster identities, too:
load clusters_n_k5stationclusters_BETTER_plus_NYC_6thcluster
load eventtimesoutput_343

%Need "dates" vector from FZRA_SynopticWeatherTyping first before running
%this!

%Get rid of first 1976-1978 records:
event_pct_lightFZRA(1:length(event_pct_lightFZRA)-length(dates)) = [];
event_stationcounts(1:length(event_stationcounts)-length(dates),:) = [];

%% Do this thing
figure(101)
plot(dates,event_pct_lightFZRA)

%Okay. Now we want it by station:

%Loop through events and record the event_pct at each station for that
%event. Effective percentage of hours that are light fzra
pct_light = NaN*zeros(size(event_stationcounts));
for event = 1:length(event_pct_lightFZRA)
    %Record the percentage of hours that are light for that event at each
    %station that registers an event:
    pct_light(event,~~event_stationcounts(event,:)) = event_pct_lightFZRA(event); 
end
    
%Now package this into a 2014-1979 yr x 97 station matrix:
pct_light_yrs = zeros(36,97);
for yearz = 1979:2014
    pct_light_yrs(yearz-1978,:) = mean(pct_light(year(dates)==yearz,:),'omitnan');
end


%And now let's average these time series into only six--one for each
%region!
for w = 1:6 %there are five clusters
    pct_light_time_series_regional_avg(w,:) = mean(pct_light_yrs(:,IDX_stations == w)','omitnan');
    
end

%OMG the time hath come! Plot 'em and hope for the best!
figure(41)
for w = 1:6
    plot(1979:2014,1-pct_light_time_series_regional_avg(w,:),'color',colorz(w,:),'LineWidth',2)
    hold on
    %Now calculate a least squares fit for the monthly time series and plot that, too:
%     coeffs(w,:) = polyfit(1979:2014,pct_light_time_series_regional_avg(w,:),1);
%     bestfit(w,:) = polyval(coeffs(w,:),1979:2014);
%     corrc = corrcoef(bestfit(w,:),pct_light_time_series_regional_avg(w,:))
%     rsquared(w) = corrc(1,2).^2;
%     
%     plot(1979:2014,bestfit(w,:),'color',colorz(w,:),'LineWidth',6)

    %We get a really steep intensification that is almost certainly due to
    %the transition to the icing sensor in 1995
    %Instead, do two different trend lines: one for pre- and post-automated
    %observations:
    coeffs_man(w,:) = polyfit(1979:1995,pct_light_time_series_regional_avg(w,1:17),1);
    bestfit_man(w,:) = polyval(coeffs_man(w,:),1979:1995);
    corrc_man = corrcoef(bestfit_man(w,:),pct_light_time_series_regional_avg(w,1:17))
    rsquared_man(w) = corrc_man(1,2).^2;
    
    %Now make trend lines for the post-sensor period:
    coeffs_auto(w,:) = polyfit(1996:2014,pct_light_time_series_regional_avg(w,18:end),1);
    bestfit_auto(w,:) = polyval(coeffs_auto(w,:),1996:2014);
    corrc_auto = corrcoef(bestfit_auto(w,:),pct_light_time_series_regional_avg(w,18:end))
    rsquared_auto(w) = corrc_auto(1,2).^2;
    
    plot(1979:1995,1-bestfit_man(w,:),'color',colorz(w,:),'LineWidth',6)
    plot(1996:2014,1-bestfit_auto(w,:),'color',colorz(w,:),'LineWidth',6)
    
    
    stats_man = fitlm(1979:1995,1-pct_light_time_series_regional_avg(w,1:17),'linear');
    stats_auto = fitlm(1996:2014,1-pct_light_time_series_regional_avg(w,18:end),'linear');
    pvals_man(w) = table2array(stats_man.Coefficients(1,4));
    pvals_auto(w) = table2array(stats_auto.Coefficients(1,4));
end

xlabel('Year')
ylabel('% FZRA Hours "Light" Intensity')
%legend('Northwest','Appalachia','Eastern Hotspots','S.-Central Canada','South Central','NYC')

    
    
statty =  fitlm(1996:2014,randn(19,1),'linear')






%% Nonevent trend analysis (for now, run FZRA_EventTimes from SynopticWeatherTyping_withTwoVariables and then come over here for the analysis.
for yr = 1976:2014
    nonevents_yearly(yr-1975) = sum(year(nonevent_times) == yr);
end

coeffs_non = polyfit(1976:2014,nonevents_yearly,1);
bestfit_non = polyval(coeffs_non,1976:2014);
corrc_non = corrcoef(bestfit_non,nonevents_yearly)
rsquared_non = corrc_non(1,2).^2;
stats_non = fitlm(1976:2014,nonevents_yearly,'linear');
pvals_non = table2array(stats_non.Coefficients(1,4))
   
    
figure(10)
plot(1976:2014,nonevents_yearly)
hold on
plot(1976:2014,bestfit_non)
xlabel('Year')
ylabel('No. Nonevent Reports Regionwide')

    
