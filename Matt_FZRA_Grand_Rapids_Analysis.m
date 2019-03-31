%FZRA Analysis of Grand Rapids for Laura

load a_mi.mat

%% Plot historical yearly avg trend at GRR:
%Make a matrix of hourly reports with rows as years and columns as hours:
GRRHourlyFreq_yearly = squeeze(a_mi.HourlyFreq_yearly(7,:,:));

t = 1976:2014;
y = sum(GRRHourlyFreq_yearly,2)';
figure(1)
plot(t,y)
xlabel('Year')
ylabel('Hourly FZRA reports at GRR')
title('Yearly Total Hours of FZRA Reported at Grand Rapids Airport')
hold on
p = polyfit(t,y,1);
yfit = polyval(p,t);
plot(t,yfit,'k:')
legend('Annual Total','Linear Fit')
grid on

[r,pval] = corrcoef(t,y);

stats = regstats(1976:2014,sum(GRRHourlyFreq_yearly,2)')


%% Plot the time-of-day trends at GRR:
figure(5)
plot(0:23,sum(HourlyFreq)/sum(sum(HourlyFreq))*100)
ylabel('Frequency (%)')
xlabel('Hour (UTC)')


%Let's compare trends at GR:
%chooseyears = 1976:1985;
%chooseyears = 1986:1995;
%chooseyears = 1996:2005;
%chooseyears = 2006:2014;
chooseyears = 1976:2014;
%hourlytrends = squeeze(sum(GRRHourlyFreq_yearly),1));
hourlytrends = sum(GRRHourlyFreq_yearly(chooseyears - 1975, :));
hourlytrends5 = hourlytrends/sum(hourlytrends)*100;

load hourlytrendsGR 

figure(3)
plot(0:23,hourlytrends5)
hold on
plot(0:23,hourlytrends1)
plot(0:23,hourlytrends2)
plot(0:23,hourlytrends3)
plot(0:23,hourlytrends4)
ylabel('Frequency (%)')
xlabel('Hour (UTC)')
title('Breakdown of FZRA Reports by Time of Day at Grand Rapids')
grid on

legend('Average','1976 - 1985','1986 - 1995','1996-2005','2006 - 2014')


%% Non-parametric Significance Analysis with Mann-Kendall
tic
[taub tau h sig Z S sigma sen n senplot cilower ciupper] = ktaub([t;y]',0.6,1);
toc

