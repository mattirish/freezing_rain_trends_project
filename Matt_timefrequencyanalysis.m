%% Matt_timefrequencyanalysis.m
% March 2016
clear 

load a_allstations.mat
load numobsperyear.mat

HourlyFreq = [a_ia.HourlyFreq;
                    a_il.HourlyFreq;
                    a_in.HourlyFreq;
                    a_mi.HourlyFreq;
                    a_mn.HourlyFreq;
                    a_ny.HourlyFreq;
                    a_oh.HourlyFreq;
                    a_onqb.HourlyFreq;
                    a_pa.HourlyFreq;
                    a_wi.HourlyFreq];
              
HourlyFreq_yearly = [a_ia.HourlyFreq_yearly;
                    a_il.HourlyFreq_yearly;
                    a_in.HourlyFreq_yearly;
                    a_mi.HourlyFreq_yearly;
                    a_mn.HourlyFreq_yearly;
                    a_ny.HourlyFreq_yearly;
                    a_oh.HourlyFreq_yearly;
                    a_onqb.HourlyFreq_yearly;
                    a_pa.HourlyFreq_yearly;
                    a_wi.HourlyFreq_yearly];

%Plot the hourly frequency of hourly reports:

figure(1)
plot(0:23,sum(HourlyFreq)/sum(sum(HourlyFreq))*100)
ylabel('Frequency (%)')
xlabel('Hour (UTC)')

%Now let's see how that hourly trend has changed over time:
trends = squeeze(sum(HourlyFreq_yearly,1));

figure(2)
ylabel('Frequency (%)')
xlabel('Hour (UTC)')
for chooseyears = 1976:2014;
    hourlytrends = trends(chooseyears - 1975, :);
    hourlytrends = hourlytrends/sum(hourlytrends)*100;
    
    plot(0:23,hourlytrends)

    pause(.2)
end

%Wow, ridiculous variability. Let's do the first period and get out of
%here:
chooseyears = 2005:2014;
hourlytrends = squeeze(sum(HourlyFreq_yearly,1));
hourlytrends = sum(hourlytrends(chooseyears - 1975, :));
hourlytrends = hourlytrends/sum(hourlytrends)*100;

figure(3)
plot(0:23,hourlytrends1)
hold on
plot(0:23,hourlytrends2)
plot(0:23,hourlytrends3)
plot(0:23,hourlytrends4)
ylabel('Frequency (%)')
xlabel('Hour (UTC)')
title('Frequency of FZRN Reports by Time of Day')
grid on

legend('1976 - 1990','1991 - 2004','2005-2014','1976 - 2014')





