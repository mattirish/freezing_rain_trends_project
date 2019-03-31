%% Matt_spatialanalysis.m
% February 2016

load a_allstations2.mat %contains all the frequency and spatial data for each region as structures
load numobsperyear.mat %contains number of observations per year at each station for each region

%Concatenate all the stations from each region to create... A MEGAMATRIX:
StationLocations = [a_ia.StationLocations;
                    a_il.StationLocations;
                    a_in.StationLocations;
                    a_mi.StationLocations;
                    a_mn.StationLocations;
                    a_ny.StationLocations;
                    a_oh.StationLocations;
                    a_onqb.StationLocations;
                    a_pa.StationLocations;
                    a_wi.StationLocations];

WindSpeed = [a_ia.WindSpeed ...
                    a_il.WindSpeed ...
                    a_in.WindSpeed ...
                    a_mi.WindSpeed ...
                    a_mn.WindSpeed ...
                    a_ny.WindSpeed ...
                    a_oh.WindSpeed ...
                    a_onqb.WindSpeed ...
                    a_pa.WindSpeed ...
                    a_wi.WindSpeed];

WindDir = [a_ia.WindDir ...
                    a_il.WindDir ...
                    a_in.WindDir ...
                    a_mi.WindDir ...
                    a_mn.WindDir ...
                    a_ny.WindDir ...
                    a_oh.WindDir ...
                    a_onqb.WindDir ...
                    a_pa.WindDir ...
                    a_wi.WindDir];

%% Initial mapping and wind map
%Specify a map area of interest (the bounds of our map):
lat_lim = [38 50]; %deg N
lon_lim = [-95 -71]; %deg W

%Import the shapefile:
states = shaperead('usastatehi','UseGeoCoords',true);
provinces = shaperead('province','UseGeoCoords',true);

%Let's plot this along with the stations to make sure we've got the coastline situated correctly.
figure(1)
worldmap(lat_lim,lon_lim) %plots empty axes
title('Wind During FZRN Events in the Great Lakes');
geoshow(states,'FaceColor',[.5 1 .5])
geoshow(provinces,'FaceColor',[.5 1 .5])
framem('ffacecolor',[.5 .7 .9]); %shows water as blue
plotm(StationLocations(:,1),StationLocations(:,2),'r+') %plots measurement stations with red star
%Fix the wind dirs:
quiverm(StationLocations(:,1)',StationLocations(:,2)',cosd(-WindDir).*WindSpeed,sind(-WindDir).*WindSpeed,0.05);
xlabel('Longitude (deg W)')
ylabel('Latitude (deg N)')
title('')

%% 2D-Interpolating between points to make a raster map of average freezing rain freq
%We can use the griddata function to interpolate the data using
%the cubic method (better than linear interpolation for modeling
%things that vary smoothly across space)
%Initialize domain variables, etc.
resolution = .2;
x = (lon_lim(1):resolution:lon_lim(2));
y = (lat_lim(1):resolution:lat_lim(2));
[domain_y, domain_x] = meshgrid(x,y);

%Normalize to the total number of hours reported at each station:
numobservations = [numobsperyear_ia;
                    numobsperyear_il;
                    numobsperyear_in;
                    numobsperyear_mi;
                    numobsperyear_mn;
                    numobsperyear_ny;
                    numobsperyear_oh;
                    numobsperyear_onqb;
                    numobsperyear_pa;
                    numobsperyear_wi];

YearFreq = [a_ia.YearFreq;
                    a_il.YearFreq;
                    a_in.YearFreq;
                    a_mi.YearFreq;
                    a_mn.YearFreq;
                    a_ny.YearFreq;
                    a_oh.YearFreq;
                    a_onqb.YearFreq;
                    a_pa.YearFreq;
                    a_wi.YearFreq];
                
MonthFreq = [a_ia.MonthFreq;
                    a_il.MonthFreq;
                    a_in.MonthFreq;
                    a_mi.MonthFreq;
                    a_mn.MonthFreq;
                    a_ny.MonthFreq;
                    a_oh.MonthFreq;
                    a_onqb.MonthFreq;
                    a_pa.MonthFreq;
                    a_wi.MonthFreq];
              
%Find relative frequencies by dividing num FZRN hours recorded each year per station by number of hours recorded                
YearFreqrel = YearFreq./numobservations;
YearFreqrel = YearFreqrel*100*100; %scale by 100, then make percent

%numobservations = frachoursavail*350264; %scale fracs by total hours btwn 1976 and end of 2014
%YearFreqrel = bsxfun(@times,YearFreq,1./(numobservations./39)); %divide each col of YearFreq by vec of est total num observations
%YearFreqrel = YearFreqrel*100*100; %scale by 100, then make percent

%Also include the four corner points of our region of interest with average
%frequencies so that the interpolation is carried out on the whole region:
stationlocations_wboundpts = [StationLocations; lat_lim' lon_lim'; lat_lim' fliplr(lon_lim)'];
MonthFreq_wboundpts = [MonthFreq; ones(1,12).*mean(MonthFreq); ones(1,12).*mean(MonthFreq); ...
    ones(1,12).*mean(MonthFreq); ones(1,12).*mean(MonthFreq)];
YearFreqrel_wboundpts = [YearFreqrel; ones(1,39).*mean(YearFreqrel); ones(1,39).*mean(YearFreqrel); ...
    ones(1,39).*mean(YearFreqrel); ones(1,39).*mean(YearFreqrel)];

%Set the boundaries to be near the local values:
MonthFreq_wboundpts = [MonthFreq; ones(1,12).*mean(MonthFreq); ones(1,12).*mean(MonthFreq); ...
    ones(1,12).*mean(MonthFreq); ones(1,12).*mean(MonthFreq)];
YearFreqrel_wboundpts = [YearFreqrel; ones(1,39).*mean(YearFreqrel()); (YearFreqrel(75,:)+YearFreqrel(73,:))/2; ...
    ones(1,39).*mean(YearFreqrel()); (YearFreqrel(46,:)+YearFreqrel(81,:))/2];

% %Let's try sorting the matrices to get griddata to work (prob useless):
% [stationlocations_wboundpts_sorted,index] = sortrows(stationlocations_wboundpts);
% YearFreqrel_wboundpts_sorted = YearFreqrel_wboundpts(index,:);


%It's time! Interpolate the scattered data points using the griddata func:
%yearstoplot = 1976:2014; %select range of years to map average of
%yearstoplot = 1976:1990;
yearstoplot = 1990:2004;
%yearstoplot = 2005:2014;

%Set NaNs to the averages of their stations and figure out why they're NaNs
%later:
for k = 1:length(YearFreqrel_wboundpts)
    YearFreqrel_wboundpts(k,isnan(YearFreqrel_wboundpts(k,:))) = mean(YearFreqrel_wboundpts(k,~isnan(YearFreqrel_wboundpts(k,:))));
end
    
fzravals_wboundpts = mean(YearFreqrel_wboundpts(:,2005:2014 - 1975),2);
fzramap = griddata(stationlocations_wboundpts(:,1),stationlocations_wboundpts(:,2), ...
    fzravals_wboundpts,domain_x,domain_y,'natural'); %beautiful! Last parameter describes the method of interpolation.

%% DIFFMAP:
yearstoplot1 = 1976:1995;
yearstoplot2 = 1996:2014;
fzravals_wboundpts1 = mean(YearFreqrel_wboundpts(:,yearstoplot1 - 1975),2);
fzravals_wboundpts2 = mean(YearFreqrel_wboundpts(:,yearstoplot2 - 1975),2);
fzramap1 = griddata(stationlocations_wboundpts(:,1),stationlocations_wboundpts(:,2), ...
    fzravals_wboundpts1,domain_x,domain_y,'natural');
fzramap2 = griddata(stationlocations_wboundpts(:,1),stationlocations_wboundpts(:,2), ...
    fzravals_wboundpts2,domain_x,domain_y,'natural');
fzramap = fzramap2 - fzramap1;
%Convert to diff in number of hourly fzra reports per year:
fzramap = fzramap*.8760;

%% Relative diffmap. Have to adjust boundary points to average change.
yearstoplot1 = 1976:1994;
yearstoplot2 = 1995:2014;
fzravals_wboundpts1 = mean(YearFreqrel_wboundpts(:,yearstoplot1 - 1975),2);
fzravals_wboundpts2 = mean(YearFreqrel_wboundpts(:,yearstoplot2 - 1975),2);
fzramap1 = griddata(stationlocations_wboundpts(:,1),stationlocations_wboundpts(:,2), ...
    fzravals_wboundpts1,domain_x,domain_y,'natural');
fzramap2 = griddata(stationlocations_wboundpts(:,1),stationlocations_wboundpts(:,2), ...
    fzravals_wboundpts2,domain_x,domain_y,'natural');
fzramap = (fzramap2 - fzramap1)./fzramap1.*100;
%Convert to diff in number of hourly fzra reports per year:
fzramap = fzramap*.8760;


%% Plot the diffmap!
%fzrnmap19902004 = griddata(stationlocations_wboundpts(:,1),stationlocations_wboundpts(:,2),fzrnvals_wboundpts,domain_x,domain_y,'natural'); %beautiful! Last parameter describes the method of interpolation.
%fzrnmapdiff = fzrnmap - fzrnmap19902004;

% %Now we need a boundary matrix to mask the map with.
% %Convert the vector boundaries of US to a matrix using the vec2mtx function:
% cellsperdegree = 1/resolution; %resolution of map raster
% [boundary_matrix, R_dummy] = vec2mtx(states(22).Lat, states(22).Lon, cellsperdegree, lat_lim, lon_lim, 'filled');
% 
% %Add one row and column to the US boundary matrix with its "out of
% %bounds" default value of 2 to make it match the size of the raster matrix
% us_boundary = 2*ones(size(fzrnmap));
% us_boundary(2:end,2:end) = boundary_matrix;

%Plot the map! She's gonna be beautiful.
figure(2)
worldmap(lat_lim,lon_lim)
colormap(parula)
caxis([-50,50])
%demcmap([min(fzrnmap),100],64,'window','window') %applies colormap to the raster data
fzra_ref = georasterref('RasterSize', size(fzramap), 'Latlim', lat_lim, 'Lonlim', lon_lim);
geoshow(states,'FaceColor',[1 1 1])
geoshow(provinces,'FaceColor',[1 1 1])
%framem('ffacecolor',[.5 .7 .9]); %shows water as blue
geoshow(fzramap,fzra_ref,'DisplayType','texturemap','FaceAlpha',0.8)
plotm(stationlocations_wboundpts(:,1),stationlocations_wboundpts(:,2),'r+') %plots measurement stations
%caxis([0, 25])
colorbar
hold on
title('Change in No. Hourly FZRA Reports per Year from 1976-1985 to 2000-2014 Period Avgs.')
%geoshow(fzrnmap,fzrn_ref,'DisplayType','contour')
%[C h] = contourm(fzrnmap,fzrn_ref,'LevelStep',5,'LineWidth',1,'LineColor','k');
%clabelm(C)


%% Plot a bubble difference map instead!
yearstoplot1 = 1976:1994;
yearstoplot2 = 1995:2014;

figure(3)
worldmap(lat_lim,lon_lim)
colormap(parula)
caxis([-40, 40])
%demcmap([-100,100],64,'window','window') %applies colormap to the raster data
geoshow(states,'FaceColor',[1 1 1])
hold on
geoshow(provinces,'FaceColor',[1 1 1])
framem('ffacecolor',[.5 .7 .9]); %shows water as blue
%geoshow(fzrnmap,fzrn_ref,'DisplayType','texturemap','FaceAlpha',0.8)
plotm(StationLocations(:,1),StationLocations(:,2),'r+') %plots measurement stations

sizecoeff = 7;
freqs1 = mean(YearFreqrel(:,yearstoplot1 - 1975),2,'omitnan');
freqs2 = mean(YearFreqrel(:,yearstoplot2 - 1975),2,'omitnan');
circleareas = (freqs2-freqs1)./freqs1*100;
scatterm(StationLocations(:,1),StationLocations(:,2),250,circleareas,'filled')
%caxis([0, 25])
colorbar
hold on
title('Difference FZRA Reports per Year from 1976-1995 to 1996-2014 Period Avgs.')

title('')

%% Plot a bubble trends map where bubbles are sized by statistical significance!
t = 2000:2014;

for m = 1:97
    p = polyfit(t,YearFreqrel(m,t-1975),1);
    hoursperyear(m) = p(1);
    [r,pval] = corrcoef(t,YearFreqrel(m,t-1975));
    pvals(m) = pval(1,2);
end


figure(4)
worldmap(lat_lim,lon_lim)
colormap(parula)
caxis([-1, 1])
geoshow(states,'FaceColor',[1 1 1])
hold on
geoshow(provinces,'FaceColor',[1 1 1])
framem('ffacecolor',[.5 .7 .9]); %shows water as blue
%geoshow(fzramap,fzra_ref,'DisplayType','texturemap','FaceAlpha',0.8)
plotm(StationLocations(:,1),StationLocations(:,2),'r+') %plots measurement stations

sizecoeff = 700;
scatterm(StationLocations(:,1),StationLocations(:,2),400*(1-pvals),hoursperyear,'filled')
%caxis([0, 25])
colorbar
hold on
title('Linear Trend in Freezing Rain (2000-2014)')


%% Statistical significance analysis:

%Plot historical trends for places where the trend was more than 1 hr/year:
figure (101)
plot(1976:2014,YearFreqrel(hoursperyear > 1,:))

%Less than 1 hr/year:
figure (102)
plot(1976:2014,YearFreqrel(hoursperyear < -1,:))

%Plot trends for several Canadian stations where very significant increases in FZRA for 2000:2014:
figure (103)
plot(1976:2014,YearFreqrel(80:81,:))

%Do a Mann-Kendall test for data in the south-central Canada region where
%we see an increase. If we do one test for several stations what happens to
%our trend and its significance?
t = 1976:2014
sigincreases = YearFreqrel([75 80 82],t-1975)
sigincreases = sigincreases(:)'
tmulti = [t t t]
[taub tau h sig Z S sigma sen n senplot cilower ciupper] = ktaub([tmulti;sigincreases]',0.3,1)


%Try the mean or median station observation for each year instead:
t = 1976:2014
sigincreases = YearFreqrel([75 80 82],t-1975)
sigincreases = mean(sigincreases)
[taub tau h sig Z S sigma sen n senplot cilower ciupper] = ktaub([t;sigincreases]',0.3,1)
%We get a p-value of 0.087 if we average the three stations. Mayb we should
%fig out a way of doing this by weighting spatially so we could arrive at a
%"hotspot" of statistical certainty between these three stations, for
%example.

%% Boxplot   
% I need to normalize these by year but I don't have annual data here
% yet--update this!

%Divide by 39 years of data and avg number of hours in a month:
avghrsinamonth = 30.4*24;
for k = 1: length(MonthFreq)
    MonthFreqrel(k,:) = MonthFreq(k,:)./(39*avghrsinamonth)/.9743; %.9743 is to correct for num observations
end

figure(10)
boxplot([MonthFreqrel(:,10:12)*100 MonthFreqrel(:,1:5)*100])
set(gca, 'XTick',1:8, 'XTickLabel',{'Oct' 'Nov' 'Dec' 'Jan' 'Feb' 'Mar' 'Apr' 'May'})
ylabel('Relative Frequency (%)')
grid on


%% Mann-Kendall Interpolation
t = 1985:2014;

for m = 1:97
    tic
    [taub(m) tau(m) h(m) sig(m) Z(m) S(m) sigma(m) sen(m) n senplot cilower(m) ciupper(m)] = ktaub([t;YearFreqrel(m,t-1975)]',0.05,1);
    toc
end


%Set the boundaries to be near the local values:
sen_wboundpts = [sen'; mean(sen); (sen(75)+sen(73))/2; ...
    mean(sen); (sen(46)+sen(81))/2];
sig_wboundpts = [sig'; mean(sig); (sig(75)+sig(73))/2; ...
    mean(sig); (sig(46)+sig(81))/2];

senmap = griddata(stationlocations_wboundpts(:,1),stationlocations_wboundpts(:,2), ...
    sen_wboundpts,domain_x,domain_y,'natural'); %beautiful! Last parameter describes the method of interpolation.

sigmap = griddata(stationlocations_wboundpts(:,1),stationlocations_wboundpts(:,2), ...
    sig_wboundpts,domain_x,domain_y,'natural'); %beautiful! Last parameter describes the method of interpolation.
sigmap(find(sigmap > 0.3)) = NaN;

figure(11)
worldmap(lat_lim,lon_lim)
colormap(parula)
%caxis([-1, 1])
sen_ref = georasterref('RasterSize', size(senmap), 'Latlim', lat_lim, 'Lonlim', lon_lim);
geoshow(states,'FaceColor',[1 1 1])
hold on
geoshow(provinces,'FaceColor',[1 1 1])
framem('ffacecolor',[.5 .7 .9]); %shows water as blue
geoshow(senmap,sen_ref,'DisplayType','texturemap','FaceAlpha',0.8)
colorbar
geoshow(sigmap,sen_ref,'DisplayType','texturemap','FaceAlpha',0.4)
plotm(StationLocations(:,1),StationLocations(:,2),'r+') %plots measurement stations

hold on
title('Linear Trend in Freezing Rain (2000-2014)')


%% SEASONAL Mann-Kendall Interpolation
t = 1976:2014;

%Time domain for all months of each year 1976 to 2014 (not v elegant):
for k = 1:length(t)
    tmonths((k+(11*(k-1))):(k+(11*(k-1)))+11) = 1:12;
    tyears((k+(11*(k-1))):(k+(11*(k-1)))+11) = (k+1975)*ones(1,12);
end

%Now we just need to make an actual monthly time series and we can compute
%the seasonal Mann-Kendall. 
[taubsea tausea Sens h sig sigAdj Zs Zmod Ss Sigmas CIlower CIupper] ...
                = sktt([tyears;tmonths;MonthFreqrel(75,t-1975)]',0.3,0,1)
                
for m = 1:97
    tic
    [taub(m) tau(m) h(m) sig(m) Z(m) S(m) sigma(m) sen(m) n senplot cilower(m) ciupper(m)] = ktaub([t;MonthFreqrel(m,t-1975)]',0.05,1);
    toc
end


%Set the boundaries to be near the local values:
sen_wboundpts = [sen'; mean(sen); (sen(75)+sen(73))/2; ...
    mean(sen); (sen(46)+sen(81))/2];
sig_wboundpts = [sig'; mean(sig); (sig(75)+sig(73))/2; ...
    mean(sig); (sig(46)+sig(81))/2];

senmap = griddata(stationlocations_wboundpts(:,1),stationlocations_wboundpts(:,2), ...
    sen_wboundpts,domain_x,domain_y,'natural'); %beautiful! Last parameter describes the method of interpolation.

sigmap = griddata(stationlocations_wboundpts(:,1),stationlocations_wboundpts(:,2), ...
    sig_wboundpts,domain_x,domain_y,'natural'); %beautiful! Last parameter describes the method of interpolation.
sigmap(find(sigmap > 0.3)) = NaN;

figure(11)
worldmap(lat_lim,lon_lim)
colormap(parula)
%caxis([-1, 1])
sen_ref = georasterref('RasterSize', size(senmap), 'Latlim', lat_lim, 'Lonlim', lon_lim);
geoshow(states,'FaceColor',[1 1 1])
hold on
geoshow(provinces,'FaceColor',[1 1 1])
framem('ffacecolor',[.5 .7 .9]); %shows water as blue
geoshow(senmap,sen_ref,'DisplayType','texturemap','FaceAlpha',0.8)
colorbar
geoshow(sigmap,sen_ref,'DisplayType','texturemap','FaceAlpha',0.4)
plotm(StationLocations(:,1),StationLocations(:,2),'r+') %plots measurement stations

hold on
title('Linear Trend in Freezing Rain (2000-2014)')
