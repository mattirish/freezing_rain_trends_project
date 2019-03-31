clear
clf
% Getting the data needed into vectors
data = ['OriginalAlpena-726390         ';'OriginalBC-725396             ';'OriginalDetroitColeman-725375 ';'OriginalDetroitWillow-725376  ';...
    'OriginalDTW-725370            ';'OriginalFlint-726370          ';'OriginalGR-726350             ';'OriginalHancock-727440        ';'OriginalHoughtonL-726380      ';'OriginalIronM-727437          ';...
    'OriginalJackson-725395        ';'OriginalKalamazoo-726357      ';'OriginalLansing-725390        ';'OriginalMuskegon-726360       ';'OriginalPellston-727347       ';...
    'OriginalPontiac-726375        ';'OriginalSaginaw-726379        ';'OriginalSaultSM-727340        ';'OriginalTraverseC-726387      '];
data = cellstr(data);
YearFreq = zeros(20,39);
MonthFreq = zeros(20,12);

for z = 1:length(data)
    DATA = char(data(z));
    num = str2num(DATA((length(DATA)-5):length(DATA)));
    data2(z,1:29) = sprintf('FreezingRainFilter_%i.txt',num);
end


for i = 1:length(data)
    
Format = {'%f%f%24c','HeaderLines',1,'Delimiter' , '\r\n'};
fid = fopen(data2(i,1:29),'r');
Data = textscan(fid, Format{:});
Year = cell2mat(Data(1));
Year2(i,1:length(Year)) = cell2mat(Data(1));
Month2(i,1:length(Year)) = cell2mat(Data(2));
Month = cell2mat(Data(2));

% Monthly Proportions? Monthly Frequencies (Month Occurence/(total hours in Month))
% Monthly Frequencies
% Month(1) corresponds to January, Month(2) to February... Year(12) to December
for j = 1:12
    % Finding total freezing rain events over full time span
    MonthFreq(i,j) = numel(find(Month == j,length(Month))); % Rows is Stations, Columns is Months (Jan-Dec)
end
% numel(Aug) == 0, Aug is empty, so can use

% Yearly Frequencies
% Year(1) corresponds to 1976, Year(2) to 1977... Year(39) to 2014

for k = 1:39
    YearFreq(i,k) = numel(find(Year == (1975 + k),length(Year))); % Rows is Stations, Columns is Years (1976-2014)
end

clear DATA Format fid Data


end

%% First Look Plots at Data, trying to get a feel for what we have
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Average Freezing Rain events per certain month
MonthFreqAvgm = mean(MonthFreq)/39;

% Average Freezing Rain events per hour (according to month)
MonthFreqAvgh(1) = MonthFreqAvgm(1)/(31*24);
MonthFreqAvgh(2) = MonthFreqAvgm(2)/((((28*29)+(29*10))/39)*24);
MonthFreqAvgh(3) = MonthFreqAvgm(3)/(31*24);
MonthFreqAvgh(4) = MonthFreqAvgm(4)/(30*24);
MonthFreqAvgh(5) = MonthFreqAvgm(5)/(31*24);
MonthFreqAvgh(6) = MonthFreqAvgm(6)/(30*24);
MonthFreqAvgh(7) = MonthFreqAvgm(7)/(31*24);
MonthFreqAvgh(8) = MonthFreqAvgm(8)/(31*24);
MonthFreqAvgh(9) = MonthFreqAvgm(9)/(30*24);
MonthFreqAvgh(10) = MonthFreqAvgm(10)/(31*24);
MonthFreqAvgh(11) = MonthFreqAvgm(11)/(30*24);
MonthFreqAvgh(12) = MonthFreqAvgm(12)/(31*24);

% Average Freezing Rain events per year 
YearFreqAvg = mean(YearFreq);

% Also need to plot the frequency per station I have in the table
% Also look into plotting individual stations
Months = cellstr(['Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec']);

figure(1)
plot(1:12,MonthFreqAvgm,'k','linewidth',2)
xlim([1,12])
set(gca,'Xtick',1:12)
set(gca,'XTickLabel',Months')
xlabel('Month')
ylabel('Freezing Rain Events / Month')
title('Average Monthly Freezing Rain Events')

figure(2)
plot(1:12,MonthFreqAvgh,'k','linewidth',2)
xlim([1,12])
set(gca,'Xtick',1:12)
set(gca,'XTickLabel',Months')
xlabel('Month')
ylabel('Freezing Rain Events / hour')
title('Average Monthly Freezing Rain Events')

% Yearly Averages over all stations
figure(3)
plot(1976:2014,YearFreqAvg,'k','linewidth',2)
xlabel('Year')
ylabel('Freezing Rain Events/Year')
title('Average Yearly Freezing Rain Events')

% Comparing low, mid, high frequency stations
figure(4)
plot(1:12,MonthFreq(13,1:12)/39,'k','linewidth',2)
hold on 
plot(1:12,MonthFreq(10,1:12)/39,'linewidth',2)
plot(1:12,MonthFreq(12,1:12)/39,'k--','linewidth',2)
xlim([1,12])
set(gca,'Xtick',1:12)
set(gca,'XTickLabel',Months')
xlabel('Month')
ylabel('Freezing Rain Events / Month')
title('Average Monthly Freezing Rain Events')
legend('Lansing','Iron Mountain','Kalamzoo')

% Yearly Totals (all stations)
totalevent = sum(YearFreq);
figure(5)
plot(1976:2014,totalevent,'k','linewidth',2)
xlabel('Year')
ylabel('Ttal Freezing Rain Events')
title('Michigan Total Yearly Freezing Rain Events')
legend('Yearly Freezing Rain')

% Individual Stations Total Events/year
% West Side
figure(6)
plot(1976:2014,YearFreq(13,1:39),'k')
hold on
plot(1976:2014,YearFreq(8,1:39),'b')
plot(1976:2014,YearFreq(20,1:39),'r')
plot(1976:2014,YearFreq(11,1:39),'g')
plot(1976:2014,YearFreq(9,1:39),'m')
legend('Kalamzoo','Grand Rapids','Traverse City','Iron Mountain','Houghton')

figure(7)
plot(1976:2014,YearFreq(6,1:39),'k')
hold on
plot(1976:2014,YearFreq(7,1:39),'b')
plot(1976:2014,YearFreq(18,1:39),'r')
plot(1976:2014,YearFreq(1,1:39),'g')
plot(1976:2014,YearFreq(19,1:39),'m')
legend('Detroit','Flint','Saginaw','Alpena','St.Marie')

%% Looking for trends in individual months
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Conclusive Notes at the bottom
jan(1:39) = 0;
feb(1:39) = 0;
mar(1:39) = 0;
apr(1:39) = 0;
may(1:39) = 0;
sep(1:39) = 0;
oct(1:39) = 0;
nov(1:39) = 0;
dec(1:39) = 0;

for i = 1:39 % (1976 - 2014)
    for x = 1:length(data) % Stations
        for j = 1:778 % Length of Year and Month matrices
            if Year2(x,j) == 1975+i
                if Month2(x,j) == 1
                    jan(i) = jan(i) + 1; % each index is a new year, should record total occurences in jan per year
                elseif Month2(x,j) == 2
                    feb(i) = feb(i) + 1;
                elseif Month2(x,j) == 3
                    mar(i) = mar(i) + 1;
                elseif Month2(x,j) == 4
                    apr(i) = apr(i) + 1;
                elseif Month2(x,j) == 5
                    may(i) = may(i) + 1;
                elseif Month2(x,j) == 9
                    sep(i) = sep(i) + 1;
                elseif Month2(x,j) == 10
                    oct(i) = oct(i) + 1;
                elseif Month2(x,j) == 11
                    nov(i) = nov(i) + 1;
                elseif Month2(x,j) == 12
                    dec(i) = dec(i) + 1;
                end
            end
        end
    end
end


x = 1976:2014;
figure(8)

subplot(3,3,1)
plot(1976:2014,jan,'r')
xlabel('Year')
ylabel('Total Events/Year')
hold on
p = polyfit(x,jan,1);
yfit = polyval(p,x);
plot(x,yfit,'k')
legend('January','Linear Fit')

subplot(3,3,2)
plot(1976:2014,feb,'r')
xlabel('Year')
ylabel('Total Events/Year')
hold on
p = polyfit(x,feb,1);
yfit = polyval(p,x);
plot(x,yfit,'k')
legend('February','Linear Fit')
subplot(3,3,3)

plot(1976:2014,mar,'r')
xlabel('Year')
ylabel('Total Events/Year')
hold on
p = polyfit(x,mar,1);
yfit = polyval(p,x);
plot(x,yfit,'k')
legend('March','Linear Fit')

subplot(3,3,4)
plot(1976:2014,apr,'r')
xlabel('Year')
ylabel('Total Events/Year')
hold on
p = polyfit(x,apr,1);
yfit = polyval(p,x);
plot(x,yfit,'k')
legend('April','Linear Fit')

subplot(3,3,5)
plot(1976:2014,may,'r')
xlabel('Year')
ylabel('Total Events/Year')
hold on
p = polyfit(x,may,1);
yfit = polyval(p,x);
plot(x,yfit,'k')
legend('May','Linear Fit')

subplot(3,3,6)
plot(1976:2014,sep,'r')
xlabel('Year')
ylabel('Total Events/Year')
hold on
p = polyfit(x,sep,1);
yfit = polyval(p,x);
plot(x,yfit,'k')
legend('September','Linear Fit')

subplot(3,3,7)
plot(1976:2014,oct,'r')
xlabel('Year')
ylabel('Total Events/Year')
hold on
p = polyfit(x,oct,1);
yfit = polyval(p,x);
plot(x,yfit,'k')
legend('October','Linear Fit')

subplot(3,3,8)
plot(1976:2014,nov,'r')
xlabel('Year')
ylabel('Total Events/Year')
hold on
p = polyfit(x,nov,1);
yfit = polyval(p,x);
plot(x,yfit,'k')
legend('November','Linear Fit')

subplot(3,3,9)
plot(1976:2014,dec,'r')
xlabel('Year')
ylabel('Total Events/Year')
hold on
p = polyfit(x,dec,1);
yfit = polyval(p,x);
plot(x,yfit)
plot(x,yfit,'k')
legend('December','Linear Fit')

% Major Anamolie in May 2013, major storm hit Northern Michigan (UP)
% Iron Mountain and Hancock account for all of the 20 May events recorded in 2013
% https://www.facebook.com/NWSMarquette/posts/622108431152385 <<-- This is a citations of the storm I found on the web

%% Looking for different trends due to changing measuring methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Yearly Totals (all stations) different time perioods
% 1976 - 1990
% totalevent = sum(YearFreq); % Total number of events across all stations
figure(9)
plot(1976:1990,totalevent(1:15),'k','linewidth',2)
xlabel('Year')
ylabel('Total Freezing Rain Events')
title('Total Yearly Freezing Rain Events 1976-1990')

% 1990 - 2014
figure(10)
plot(1991:2014,totalevent(16:39),'k','linewidth',2)
xlabel('Year')
ylabel('Ttal Freezing Rain Events')
title('Total Yearly Freezing Rain Events 1991-2014')

% getting rid of high maximums to evaluate trend from 1991 - 2014
edittotal = totalevent;
edittotal(38) = mean(totalevent(16:39));
edittotal(32) = mean(totalevent(16:39));

figure(11)
plot(1991:2014,edittotal(16:39),'k','linewidth',2)
xlabel('Year')
ylabel('Total Freezing Rain Events')
title('Edited Total Yearly Freezing Rain Events 1991-2014')

% Doing same for the full time series
edittotal(38) = mean(totalevent(1:39));
edittotal(32) = mean(totalevent(1:39));

figure(12)
plot(1976:2014,edittotal(1:39),'k','linewidth',2)
xlabel('Year')
ylabel('Total Freezing Rain Events')
title('Edited Total Yearly Freezing Rain Events 1976-2014')
hold on 
p = polyfit(1980:2008,edittotal(5:33),1);
yfit = polyval(p,1980:2008);
plot(1980:2008,yfit,'r')
legend('Edited Yearly Totals','Linear Fit')

% 2005 - 2014
figure(13)
plot(2005:2014,totalevent(30:39),'k','linewidth',2)
xlabel('Year')
ylabel('Total Freezing Rain Events')
title('Total Yearly Freezing Rain Events 2005-2014')

% Looking at the affect of Automated Measurements
edittotal(32) = mean(totalevent(30:39));
edittotal(38) = mean(totalevent(30:39));
automateav = mean(edittotal(30:39));
automatediff = totalevent - automateav;
figure(30)
bar(1976:2014, automatediff)
xlabel('Year')
ylabel('Difference in Freezing Rain Events')
title('Michigan Impact of Automated Measurements')

bar1(1:15) = mean(totalevent(1:15)); % Mean from 1976 - 1990
bar2(1:14) = mean(totalevent(16:29)); % Mean from 1991 - 2004
bar3(1:10) = mean(totalevent(30:39)); % Mean from 2005 - 2014
bar4(1:10) = automateav; % Edited mean from 2005 - 2014 (big peaks reduced)
figure(31)
plot(1976:2014,totalevent,'k','linewidth',2)
xlabel('Year')
ylabel('Total Freezing Rain Events')
title('Michigan Total Yearly Freezing Rain Events')
hold on 
plot(1976:1990,bar1,'g','linewidth',2)
plot(1991:2004,bar2,'r','linewidth',2)
plot(2005:2014,bar3,'b','linewidth',2)
plot(2005:2014,bar4,'b--','linewidth',2)
legend('Total Events','1976-1990 Average','1991-2004 Average','2005-2014 Average','Edited 2005-2014 Average')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Separating stations by latitude in order to look at March Trend

% High Latitude (Iron Mountain, Hancock, Pellston, Alpena, Sault St. Marie)
high_jan(1:39) = 0;
high_feb(1:39) = 0;
high_mar(1:39) = 0;
high_apr(1:39) = 0;
high_may(1:39) = 0;
high_sep(1:39) = 0;
high_oct(1:39) = 0;
high_nov(1:39) = 0;
high_dec(1:39) = 0;
location = [10,8,1,18,15];

for i = 1:39 % (1976 - 2014)
    for x = 1:5 % Locations
        for j = 1:778 % Length of Year and Month matrices
            if Year2(location(x),j) == 1975+i
                if Month2(location(x),j) == 1
                    high_jan(i) = high_jan(i) + 1; % each index is a new year, should record total occurences in jan per year
                elseif Month2(location(x),j) == 2
                    high_feb(i) = high_feb(i) + 1;
                elseif Month2(location(x),j) == 3
                    high_mar(i) = high_mar(i) + 1;
                elseif Month2(location(x),j) == 4
                    high_apr(i) = high_apr(i) + 1;
                elseif Month2(location(x),j) == 5
                    high_may(i) = high_may(i) + 1;
                elseif Month2(location(x),j) == 9
                    high_sep(i) = high_sep(i) + 1;
                elseif Month2(location(x),j) == 10
                    high_oct(i) = high_oct(i) + 1;
                elseif Month2(location(x),j) == 11
                    high_nov(i) = high_nov(i) + 1;
                elseif Month2(location(x),j) == 12
                    high_dec(i) = high_dec(i) + 1;
                end
            end
        end
    end
end

% Mid Latitude (Saginaw, Flint, Grand Rapids, Muskegeon, Houghton Lake)
mid_jan(1:39) = 0;
mid_feb(1:39) = 0;
mid_mar(1:39) = 0;
mid_apr(1:39) = 0;
mid_may(1:39) = 0;
mid_sep(1:39) = 0;
mid_oct(1:39) = 0;
mid_nov(1:39) = 0;
mid_dec(1:39) = 0;
location = [17,6,7,14,9];

for i = 1:39 % (1976 - 2014)
    for x = 1:5 % Locations
        for j = 1:778 % Length of Year and Month matrices
            if Year2(location(x),j) == 1975+i
                if Month2(location(x),j) == 1
                    mid_jan(i) = mid_jan(i) + 1; % each index is a new year, should record total occurences in jan per year
                elseif Month2(location(x),j) == 2
                    mid_feb(i) = mid_feb(i) + 1;
                elseif Month2(location(x),j) == 3
                    mid_mar(i) = mid_mar(i) + 1;
                elseif Month2(location(x),j) == 4
                    mid_apr(i) = mid_apr(i) + 1;
                elseif Month2(location(x),j) == 5
                    mid_may(i) = mid_may(i) + 1;
                elseif Month2(location(x),j) == 9
                    mid_sep(i) = mid_sep(i) + 1;
                elseif Month2(location(x),j) == 10
                    mid_oct(i) = mid_oct(i) + 1;
                elseif Month2(location(x),j) == 11
                    mid_nov(i) = mid_nov(i) + 1;
                elseif Month2(location(x),j) == 12
                    mid_dec(i) = mid_dec(i) + 1;
                end
            end
        end
    end
end

% Low Latitude (DTW, Jackson, Battle Creek, Kalamazoo, Detroit Coleman, Detroit Willow)
low_jan(1:39) = 0;
low_feb(1:39) = 0;
low_mar(1:39) = 0;
low_apr(1:39) = 0;
low_may(1:39) = 0;
low_sep(1:39) = 0;
low_oct(1:39) = 0;
low_nov(1:39) = 0;
low_dec(1:39) = 0;
location = [3,4,5,11,12,2];

for i = 1:39 % (1976 - 2014)
    for x = 1:5 % Locations
        for j = 1:778 % Length of Year and Month matrices
            if Year2(location(x),j) == 1975+i
                if Month2(location(x),j) == 1
                    low_jan(i) = low_jan(i) + 1; % each index is a new year, should record total occurences in jan per year
                elseif Month2(location(x),j) == 2
                    low_feb(i) = low_feb(i) + 1;
                elseif Month2(location(x),j) == 3
                    low_mar(i) = low_mar(i) + 1;
                elseif Month2(location(x),j) == 4
                    low_apr(i) = low_apr(i) + 1;
                elseif Month2(location(x),j) == 5
                    low_may(i) = low_may(i) + 1;
                elseif Month2(location(x),j) == 9
                    low_sep(i) = low_sep(i) + 1;
                elseif Month2(location(x),j) == 10
                    low_oct(i) = low_oct(i) + 1;
                elseif Month2(location(x),j) == 11
                    low_nov(i) = low_nov(i) + 1;
                elseif Month2(location(x),j) == 12
                    low_dec(i) = low_dec(i) + 1;
                end
            end
        end
    end
end

%% Plotting the Latitude Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Low Latitude

figure(14)
x = 1976:2014;

subplot(3,3,1)
plot(1976:2014,low_jan,'r')
xlabel('Year')
ylabel('Total Events/Year')
hold on
p = polyfit(x,low_jan,1);
yfit = polyval(p,x);
plot(x,yfit,'k')
legend('January','Linear Fit')

subplot(3,3,2)
plot(1976:2014,low_feb,'r')
xlabel('Year')
ylabel('Total Events/Year')
title('Low Latitude Stations','fontweight','bold','fontsize',18);
hold on
p = polyfit(x,low_feb,1);
yfit = polyval(p,x);
plot(x,yfit,'k')
legend('February','Linear Fit')
subplot(3,3,3)

plot(1976:2014,low_mar,'r')
xlabel('Year')
ylabel('Total Events/Year')
hold on
p = polyfit(x,low_mar,1);
yfit = polyval(p,x);
plot(x,yfit,'k')
legend('March','Linear Fit')

subplot(3,3,4)
plot(1976:2014,low_apr,'r')
xlabel('Year')
ylabel('Total Events/Year')
hold on
p = polyfit(x,low_apr,1);
yfit = polyval(p,x);
plot(x,yfit,'k')
legend('April','Linear Fit')

subplot(3,3,5)
plot(1976:2014,low_may,'r')
xlabel('Year')
ylabel('Total Events/Year')
hold on
p = polyfit(x,low_may,1);
yfit = polyval(p,x);
plot(x,yfit,'k')
legend('May','Linear Fit')

subplot(3,3,6)
plot(1976:2014,low_sep,'r')
xlabel('Year')
ylabel('Total Events/Year')
hold on
p = polyfit(x,low_sep,1);
yfit = polyval(p,x);
plot(x,yfit,'k')
legend('September','Linear Fit')

subplot(3,3,7)
plot(1976:2014,low_oct,'r')
xlabel('Year')
ylabel('Total Events/Year')
hold on
p = polyfit(x,low_oct,1);
yfit = polyval(p,x);
plot(x,yfit,'k')
legend('October','Linear Fit')

subplot(3,3,8)
plot(1976:2014,low_nov,'r')
xlabel('Year')
ylabel('Total Events/Year')
hold on
p = polyfit(x,low_nov,1);
yfit = polyval(p,x);
plot(x,yfit,'k')
legend('November','Linear Fit')

subplot(3,3,9)
plot(1976:2014,low_dec,'r')
xlabel('Year')
ylabel('Total Events/Year')
hold on
p = polyfit(x,low_dec,1);
yfit = polyval(p,x);
plot(x,yfit)
plot(x,yfit,'k')
legend('December','Linear Fit')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Mid Latitude

figure(15)
x = 1976:2014;

subplot(3,3,1)
plot(1976:2014,mid_jan,'r')
xlabel('Year')
ylabel('Total Events/Year')
hold on
p = polyfit(x,mid_jan,1);
yfit = polyval(p,x);
plot(x,yfit,'k')
legend('January','Linear Fit')

subplot(3,3,2)
plot(1976:2014,mid_feb,'r')
xlabel('Year')
ylabel('Total Events/Year')
title('Mid Latitude Stations','fontweight','bold','fontsize',18);
hold on
p = polyfit(x,mid_feb,1);
yfit = polyval(p,x);
plot(x,yfit,'k')
legend('February','Linear Fit')
subplot(3,3,3)

plot(1976:2014,mid_mar,'r')
xlabel('Year')
ylabel('Total Events/Year')
hold on
p = polyfit(x,mid_mar,1);
yfit = polyval(p,x);
plot(x,yfit,'k')
legend('March','Linear Fit')

subplot(3,3,4)
plot(1976:2014,mid_apr,'r')
xlabel('Year')
ylabel('Total Events/Year')
hold on
p = polyfit(x,mid_apr,1);
yfit = polyval(p,x);
plot(x,yfit,'k')
legend('April','Linear Fit')

subplot(3,3,5)
plot(1976:2014,mid_may,'r')
xlabel('Year')
ylabel('Total Events/Year')
hold on
p = polyfit(x,mid_may,1);
yfit = polyval(p,x);
plot(x,yfit,'k')
legend('May','Linear Fit')

subplot(3,3,6)
plot(1976:2014,mid_sep,'r')
xlabel('Year')
ylabel('Total Events/Year')
hold on
p = polyfit(x,mid_sep,1);
yfit = polyval(p,x);
plot(x,yfit,'k')
legend('September','Linear Fit')

subplot(3,3,7)
plot(1976:2014,mid_oct,'r')
xlabel('Year')
ylabel('Total Events/Year')
hold on
p = polyfit(x,mid_oct,1);
yfit = polyval(p,x);
plot(x,yfit,'k')
legend('October','Linear Fit')

subplot(3,3,8)
plot(1976:2014,mid_nov,'r')
xlabel('Year')
ylabel('Total Events/Year')
hold on
p = polyfit(x,mid_nov,1);
yfit = polyval(p,x);
plot(x,yfit,'k')
legend('November','Linear Fit')

subplot(3,3,9)
plot(1976:2014,mid_dec,'r')
xlabel('Year')
ylabel('Total Events/Year')
hold on
p = polyfit(x,mid_dec,1);
yfit = polyval(p,x);
plot(x,yfit)
plot(x,yfit,'k')
legend('December','Linear Fit')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% High Latitude

figure(16)
x = 1976:2014;

subplot(3,3,1)
plot(1976:2014,high_jan,'r')
xlabel('Year')
ylabel('Total Events/Year')
hold on
p = polyfit(x,high_jan,1);
yfit = polyval(p,x);
plot(x,yfit,'k')
legend('January','Linear Fit')

subplot(3,3,2)
plot(1976:2014,high_feb,'r')
xlabel('Year')
ylabel('Total Events/Year')
title('High Latitude Stations','fontweight','bold','fontsize',18);
hold on
p = polyfit(x,high_feb,1);
yfit = polyval(p,x);
plot(x,yfit,'k')
legend('February','Linear Fit')
subplot(3,3,3)

plot(1976:2014,high_mar,'r')
xlabel('Year')
ylabel('Total Events/Year')
hold on
p = polyfit(x,high_mar,1);
yfit = polyval(p,x);
plot(x,yfit,'k')
legend('March','Linear Fit')

subplot(3,3,4)
plot(1976:2014,high_apr,'r')
xlabel('Year')
ylabel('Total Events/Year')
hold on
p = polyfit(x,high_apr,1);
yfit = polyval(p,x);
plot(x,yfit,'k')
legend('April','Linear Fit')

subplot(3,3,5)
plot(1976:2014,high_may,'r')
xlabel('Year')
ylabel('Total Events/Year')
hold on
p = polyfit(x,high_may,1);
yfit = polyval(p,x);
plot(x,yfit,'k')
legend('May','Linear Fit')

subplot(3,3,6)
plot(1976:2014,high_sep,'r')
xlabel('Year')
ylabel('Total Events/Year')
hold on
p = polyfit(x,high_sep,1);
yfit = polyval(p,x);
plot(x,yfit,'k')
legend('September','Linear Fit')

subplot(3,3,7)
plot(1976:2014,high_oct,'r')
xlabel('Year')
ylabel('Total Events/Year')
hold on
p = polyfit(x,high_oct,1);
yfit = polyval(p,x);
plot(x,yfit,'k')
legend('October','Linear Fit')

subplot(3,3,8)
plot(1976:2014,high_nov,'r')
xlabel('Year')
ylabel('Total Events/Year')
hold on
p = polyfit(x,high_nov,1);
yfit = polyval(p,x);
plot(x,yfit,'k')
legend('November','Linear Fit')

subplot(3,3,9)
plot(1976:2014,high_dec,'r')
xlabel('Year')
ylabel('Total Events/Year')
hold on
p = polyfit(x,high_dec,1);
yfit = polyval(p,x);
plot(x,yfit)
plot(x,yfit,'k')
legend('December','Linear Fit')






















