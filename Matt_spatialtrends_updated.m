%% Matt_spatialtrends_updated.m
% Updated May 2019

%For AMS journal submission:
set(0,'DefaultAxesFontName', 'Helvetica');

%load a_allstations2.mat %contains all the frequency and spatial data for each region as structures
load a_all.mat %contains a, the batch summary data for all 98 stations
%load b_fixedcoords.mat %contains b, the full dataset of FZRA events for each station.
load b.mat %contains b, the full dataset of FZRA events for each station. *fixed b file. coords are fixed too. can prob delete b_fixedcoords

load yearly_trends.mat %contains trends calculated in R using modified Mann-Kendall test


%load numobsperyear.mat %contains number of observations per year at each station for each region

%Concatenate all the stations from each region to create... A MEGAMATRIX:
% StationLocations = [a_ia.StationLocations;
%                     a_il.StationLocations;
%                     a_in.StationLocations;
%                     a_mi.StationLocations;
%                     a_mn.StationLocations;
%                     a_ny.StationLocations;
%                     a_oh.StationLocations;
%                     a_onqb.StationLocations;
%                     a_pa.StationLocations;
%                     a_wi.StationLocations];
% 
% WindSpeed = [a_ia.WindSpeed ...
%                     a_il.WindSpeed ...
%                     a_in.WindSpeed ...
%                     a_mi.WindSpeed ...
%                     a_mn.WindSpeed ...
%                     a_ny.WindSpeed ...
%                     a_oh.WindSpeed ...
%                     a_onqb.WindSpeed ...
%                     a_pa.WindSpeed ...
%                     a_wi.WindSpeed];
% 
% WindDir = [a_ia.WindDir ...
%                     a_il.WindDir ...
%                     a_in.WindDir ...
%                     a_mi.WindDir ...
%                     a_mn.WindDir ...
%                     a_ny.WindDir ...
%                     a_oh.WindDir ...
%                     a_onqb.WindDir ...
%                     a_pa.WindDir ...
%                     a_wi.WindDir];

%cell2mat(struct2cell(structfun (@(x) x.coords, b, 'UniformOutput', false)))


%% Initial mapping and wind map
%Specify a map area of interest (the bounds of our map):
lat_lim = [37.9 50]; %deg N
lon_lim = [-97 -72]; %deg W

%Import the shapefiles:
states = shaperead('usastatehi','UseGeoCoords',true);
provinces = shaperead('shapefiles/province','UseGeoCoords',true);
mexstates = shaperead('shapefiles/mexstates','UseGeoCoords',true);
lakes = shaperead('shapefiles/ne_10m_lakes_north_america','UseGeoCoords',true);
greatlakes = shaperead('shapefiles/Great_Lakes','UseGeoCoords',true);
[Z_elev, refvec_elev] = etopo(1, lat_lim, lon_lim);


%Let's plot this along with the stations to make sure we've got the coastline situated correctly.
figure(1)
ax = worldmap(lat_lim,lon_lim) %plots empty axes
title('Wind During FZRA Events in the Great Lakes');
geoshow(states,'FaceColor',[.5 1 .5])
geoshow(provinces,'FaceColor',[.5 1 .5])
framem('ffacecolor',[.5 .7 .9]); %shows water as blue
plotm(a.StationLocations(:,1),a.StationLocations(:,2),'r+') %plots measurement stations with red star
% mstruct = gcm; %can check current projection with this. worldmap makes it
% an 'eqdconic'. Maybe change out for Albers Great Lakes
%Fix the wind dirs:
%quiverm(a.StationLocations(:,1)',a.StationLocations(:,2)',cosd(-a.WindDir).*a.WindSpeed,sind(-a.WindDir).*a.WindSpeed,0.05);

%Plot the station numbers as I have them on the map:
for m = 1:length(a.StationLocations)
    textm(a.StationLocations(m,1),a.StationLocations(m,2),num2str(m))
end

xlabel('Longitude (deg W)')
ylabel('Latitude (deg N)')
title('')

%% 2D-Interpolating between points to make a raster map of average freezing rain freq
%We can use the griddata function to interpolate the data using
%the cubic method (better than linear interpolation for modeling
%things that vary smoothly across space)
%Initialize domain variables, etc.
resolution = .05;
x = (lon_lim(1):resolution:lon_lim(2));
y = (lat_lim(1):resolution:lat_lim(2));
[domain_y, domain_x] = meshgrid(x,y);

%Normalize to the total number of hours reported at each station:
% numobservations = [numobsperyear_ia;
%                     numobsperyear_il;
%                     numobsperyear_in;
%                     numobsperyear_mi;
%                     numobsperyear_mn;
%                     numobsperyear_ny;
%                     numobsperyear_oh;
%                     numobsperyear_onqb;
%                     numobsperyear_pa;
%                     numobsperyear_wi];
% 
% YearFreq = [a_ia.YearFreq;
%                     a_il.YearFreq;
%                     a_in.YearFreq;
%                     a_mi.YearFreq;
%                     a_mn.YearFreq;
%                     a_ny.YearFreq;
%                     a_oh.YearFreq;
%                     a_onqb.YearFreq;
%                     a_pa.YearFreq;
%                     a_wi.YearFreq];
%                 
% MonthFreq = [a_ia.MonthFreq;
%                     a_il.MonthFreq;
%                     a_in.MonthFreq;
%                     a_mi.MonthFreq;
%                     a_mn.MonthFreq;
%                     a_ny.MonthFreq;
%                     a_oh.MonthFreq;
%                     a_onqb.MonthFreq;
%                     a_pa.MonthFreq;
%                     a_wi.MonthFreq];
              
%Find relative frequencies by dividing num FZRA hours recorded each year per station by number of hours recorded                
%YearFreqrel = YearFreq./numobservations;
%YearFreqrel = YearFreqrel*100*100; %scale by 100, then make percent
YearFreqrel = a.YearFreq_rel; %in total hours of FZRA per year, normalized for num. observations.

%Interpolate NaNs to linear values:
for m = 1:97
    YearFreqrel(m,isnan(YearFreqrel(m,:))) = interp1(find(~isnan(YearFreqrel(m,:))), YearFreqrel(m,~isnan(YearFreqrel(m,:))), find(isnan(YearFreqrel(m,:))),'linear');
end

%YearFreqrel = YearFreqrel*100*100; %scale by 100, then make percent

%Also include the four corner points of our region of interest with average
%frequencies so that the interpolation is carried out on the whole region:
stationlocations_wboundpts = [a.StationLocations; lat_lim' lon_lim'; lat_lim' fliplr(lon_lim)'];
%MonthFreq_wboundpts = [MonthFreq; ones(1,12).*mean(MonthFreq); ones(1,12).*mean(MonthFreq); ...
%    ones(1,12).*mean(MonthFreq); ones(1,12).*mean(MonthFreq)];
YearFreqrel_wboundpts = [YearFreqrel; ones(1,39).*mean(YearFreqrel); ones(1,39).*mean(YearFreqrel); ...
    ones(1,39).*mean(YearFreqrel); ones(1,39).*mean(YearFreqrel)];

%Set the boundaries to be near the local values:
YearFreqrel_wboundpts = [YearFreqrel; ones(1,39).*mean(YearFreqrel()); (YearFreqrel(75,:)+YearFreqrel(73,:))/2; ...
    ones(1,39).*mean(YearFreqrel()); (YearFreqrel(46,:)+YearFreqrel(81,:))/2];


%It's time! Interpolate the scattered data points using the griddata func:
yearstoplot = 1976:2014; %select range of years to map average of
%yearstoplot = 1976:1990;
%yearstoplot = 1990:2004;
%yearstoplot = 2005:2014;

%Set NaNs to the averages of their stations and figure out why they're NaNs
%later:
for k = 1:length(YearFreqrel_wboundpts)
    YearFreqrel_wboundpts(k,isnan(YearFreqrel_wboundpts(k,:))) = mean(YearFreqrel_wboundpts(k,~isnan(YearFreqrel_wboundpts(k,:))));
end
    
fzravals_wboundpts = mean(YearFreqrel_wboundpts(:,yearstoplot - 1975),2);

fzramap = griddata(stationlocations_wboundpts(:,1),stationlocations_wboundpts(:,2), ...
    fzravals_wboundpts,domain_x,domain_y,'natural'); %beautiful! Last parameter describes the method of interpolation.


%% DIFFMAP:
yearstoplot1 = 1976:1995;
yearstoplot2 = 1996:2014;
yearstoplot1 = 2005:2009;
yearstoplot2 = 2009:2014;
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
%fzramap = fzramap*.8760;


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

lat_plot = [38 50]; %deg N
lon_plot = [-96 -72]; %deg W

%Plot the map! She's gonna be beautiful.
figure(2)
worldmap(lat_plot,lon_plot)
colormap(parula)
%caxis([0,25])

%cptcmap('SVS_tempanomaly', 'mapping', 'scaled');
%   colorbar; 
%demcmap([min(fzramap),100],64,'window','window') %applies colormap to the raster data
fzra_ref = georasterref('RasterSize', size(fzramap), 'Latlim', lat_lim, 'Lonlim', lon_lim);
geoshow(states,'FaceColor',[1 1 1])
geoshow(provinces,'FaceColor',[1 1 1])
%framem('ffacecolor',[.5 .7 .9]); %shows water as blue
geoshow(fzramap,fzra_ref,'DisplayType','texturemap','FaceAlpha',0.8)
plotm(stationlocations_wboundpts(:,1),stationlocations_wboundpts(:,2),'r+') %plots measurement stations
%caxis([0, 25])
geoshow(states,'FaceAlpha',0)
geoshow(provinces,'FaceAlpha',0)
c = colorbar
c.FontSize = 22
%c.FontName = 'Gill Sans MT'
hold on
%title('Change in No. Hourly FZRA Reports per Year from 1976-1985 to 2000-2014 Period Avgs.')
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

%Plot an elevation basemap for reference:
[C h] = contourm(fzrnmap,fzrn_ref,'LevelStep',5,'LineWidth',1,'LineColor','k');


plotm(a.StationLocations(:,1),a.StationLocations(:,2),'ko') %plots measurement stations
sizecoeff = 7;
freqs1 = mean(YearFreqrel(:,yearstoplot1 - 1975),2,'omitnan');
freqs2 = mean(YearFreqrel(:,yearstoplot2 - 1975),2,'omitnan');
circleareas = (freqs2-freqs1)./freqs1*100;
scatterm(a.StationLocations(:,1),a.StationLocations(:,2),250,circleareas,'filled')

%Plot sizing stations by how many storms they've participated in:
% scatterm(a.StationLocations(:,1),a.StationLocations(:,2),250,sum(event_stationcounts),'filled')

%Or size stations by how many nonevent hours they've hosted:
for m = 1:97
    total_nonevents(m) = sum(nonevent_stations == m);
end
scatterm(a.StationLocations(:,1),a.StationLocations(:,2),250,total_nonevents./sum(event_stationcounts),'filled')
caxis([0, 0.3])
colorbar

    
%caxis([0, 25])
colorbar
hold on
title('Difference FZRA Reports per Year from 1976-1995 to 1996-2014 Period Avgs.')

%caxis([0, 25])
colorbar
hold on
title('Difference FZRA Reports per Year from 1976-1995 to 1996-2014 Period Avgs.')

%% Plot a bubble trends map where bubbles are sized by statistical significance!
t = 1976:2014;

for m = 1:97
    p = polyfit(t,YearFreqrel(m,t-1975),1);
    hoursperyear(m) = p(1);
    [r,pval] = corrcoef(t,YearFreqrel(m,t-1975));
    pvals(m) = pval(1,2);
end

%Look at specific months:
% Monthchoice = 11;
% for m = 1:97
%     p = polyfit(t,a.MonthFreq_yearly_rel(m,t-1975,Monthchoice),1);
%     hoursperyear(m) = p(1);
%     [r,pval] = corrcoef(t,a.MonthFreq_yearly_rel(m,t-1975,Monthchoice));
%     pvals(m) = pval(1,2);
% end

%Ultimate freezing rain map:
figure(7)
ax1 = worldmap(lat_lim,lon_lim)
colormap(ax1, flipud(gray(256)))
caxis([-20 1400])
%Plot DEM:
elev_ax = geoshow(ax1,Z_elev, refvec_elev, 'DisplayType', 'texturemap');
geoshow(lakes,'FaceColor',[1 1 1])
geoshow(greatlakes,'FaceColor',[1 1 1])
ax2 = axes;
axis(ax2,'off')
axes(ax2)
bubb_ax = worldmap(lat_lim,lon_lim)
cb_bubb = colorbar(bubb_ax,'Location','east')
cb_elev = colorbar(ax1,'Location','west')
cb_bubb.FontSize = 12
cb_elev.FontSize = 12
xlabel(cb_elev,'Elevation (m)')
xlabel(cb_bubb,'Decadal Linear Trend (hours year^{-1} decade^{-1})')
hold on
%geoshow(fzramap,fzra_ref,'DisplayType','texturemap','FaceAlpha',0.8)
%c_bubbs = scatterm(a.StationLocations(:,1),a.StationLocations(:,2),500*(1-pvals),hoursperyear,'filled')


pvals_mmkh = yearly_trends_data(2,:);
pvals_mmkh(pvals_mmkh == 1 | isnan(pvals_mmkh)) = 0.99;
%c_bubbs = scatterm(a.StationLocations(:,1),a.StationLocations(:,2),500*(1-pvals_mmkh),yearly_trends_data(7,:)'*39./(median(YearFreqrel,2,'omitnan') - yearly_trends_data(7,:)'*19)*100,'filled') % percent change '76 to '14

c_bubbs = scatterm(a.StationLocations(:,1),a.StationLocations(:,2),500*(1-pvals_mmkh),yearly_trends_data(7,:)*10,'filled') %1lag hamed

%c_bubbs = scatterm(a.StationLocations(:,1),a.StationLocations(:,2),300,yearly_trends_3lag_data(7,:),'filled') %3lag hamed
%c_bubbs = scatterm(a.StationLocations(:,1),a.StationLocations(:,2),300,median(YearFreqrel(:,(1976:1985) - 1975),2,'omitnan'),'filled') %1975-1986 average

% % Map the 1976 and 2014 "representative" years:
% % 1976 (median at 1995 minus (19 years * slope):
% c_bubbs = scatterm(a.StationLocations(:,1),a.StationLocations(:,2),300,median(YearFreqrel,2,'omitnan') - yearly_trends_data(7,:)'*19,'filled')
% % 2014 (median at 1995 plus (19 years * slope):
% %c_bubbs = scatterm(a.StationLocations(:,1),a.StationLocations(:,2),300,median(YearFreqrel,2,'omitnan') + yearly_trends_data(7,:)'*19,'filled')
% caxis(ax2,[0 40])
% xlabel(cb_bubb,'Hours of Freezing Rain Per Year)')

%c_bubbs = scatterm(a.StationLocations(:,1),a.StationLocations(:,2),100,pvals,'filled')
%scatterm(a.StationLocations(:,1),a.StationLocations(:,2),500*(1-sig),sen_percent,'filled')
cptcmap('SVS_tempanomaly', bubb_ax, 'mapping', 'scaled');
caxis(ax2,[-max(abs(yearly_trends_data(7,:)))*10 max(abs(yearly_trends_data(7,:)))*10])
%caxis(ax2,[-max(abs(yearly_trends_3lag_data(7,:))) max(abs(yearly_trends_3lag_data(7,:)))])
caxis(ax2,[0 max(abs(median(YearFreqrel(:,(1976:1985) - 1975),2,'omitnan')))])

% For percentage change plot:
cptcmap('SVS_tempanomaly', bubb_ax, 'mapping', 'scaled');
caxis(ax2,[-200 200])
ypercentages = get(cb_bubb,'YTickLabel');
perc = repmat('%',size(ypercentages,1),1);
ypercentages = strcat(ypercentages, perc);
set(cb_bubb,'YTickLabel',ypercentages);
xlabel(cb_bubb,'Change 1976?2014')


hold on
transparentshapes = makesymbolspec('Polygon',{'Default','FaceAlpha',0});
%Plot state and province outlines:
geoshow(states,'SymbolSpec',transparentshapes)
geoshow(provinces,'SymbolSpec',transparentshapes)

plotm(a.StationLocations(:,1),a.StationLocations(:,2),'k.','MarkerSize',22)
%Add a small X to locations where p<0.05:
plotm(a.StationLocations(yearly_trends_data(2,:)<0.05,1),a.StationLocations(yearly_trends_data(2,:)<0.05,2),'kx','MarkerSize',20)




figure(1000)
plot(yearly_trends_data(2,:))



%% Test for autocorrelation for annual data. This will help determine validity of using Mann-Kendall test.
t = 1976:2014;

%Use Spearman's correlation coefficient, since it's nonparametric:
for m = 1:97
    [r,pval] = corr(YearFreqrel(m,1:end-1)',YearFreqrel(m,2:end)','type','Spearman','rows','pairwise');
    rautocorr(m) = r;
    pvals(m) = pval;
end

autocorrelated = median(rautocorr,'omitnan')

sizecoeff = 700;
scatterm(a.StationLocations(:,1),a.StationLocations(:,2),500*(1-pvals),rautocorr,'filled')


%% Statistical significance analysis:

%Plot historical trends for places where the trend was more than .2 hr/year:
figure (101)
plot(1976:2014,YearFreqrel(hoursperyear > 0.2,:))

%Less than -0.2 hr/year:
figure (102)
plot(1976:2014,YearFreqrel(hoursperyear < -0.2,:))

%Plot trends for several Canadian stations where very significant increases in FZRA for 2000:2014:
figure (103)
plot(1976:2014,YearFreqrel(80:81,:))

%Do a Mann-Kendall test for data in the south-central Canada region where
%we see an increase. If we do one test for several stations what happens to
%our trend and its significance?
t = 1976:2014
sigincreases = YearFreqrel([10 12 6],t-1975)
sigincreases = sigincreases(:)'
tmulti = [t t t]
[taub tau h sig Z S sigma sen n senplot cilower ciupper] = ktaub([tmulti;sigincreases]',0.7,1)
%[taubsea tausea Sens h sig sigAdj Zs Zmod Ss Sigmas CIlower CIupper] = sktt(datain,alpha,wantplot,StartSeason)
grid on
% _________________________________________________________________________________________
%Try the mean or median station observation for each year instead:
t = 1976:2014
sigincreases = YearFreqrel([10 12 6 5],t-1975)
sigincreases = mean(sigincreases)
[taub tau h sig Z S sigma sen n senplot cilower ciupper] = ktaub([t;sigincreases]',0.05,1)
%We get a p-value of 0.087 if we average the three stations. Mayb we should
%fig out a way of doing this by weighting spatially so we could arrive at a
%"hotspot" of statistical certainty between these three stations, for
%example.

hold on
scatter(t,YearFreqrel(10,t-1975))
scatter(t,YearFreqrel(12,t-1975))
scatter(t,YearFreqrel(6,t-1975))
scatter(t,YearFreqrel(5,t-1975))
ylabel('Total Yearly Hrs. of FZRA')


%Let's do that for the decreases in Appalachia, too:
% _________________________________________________________________________________________
%Try the mean or median station observation for each year instead:
t = 1976:2014
indices = find(a.StationLocations(:,1) > 40.1 & a.StationLocations(:,1) < 43.5 & ...
               a.StationLocations(:,2) < -75 & a.StationLocations(:,2) > -80) 
sigincreases = YearFreqrel(indices,t-1975)
sigincreases = mean(sigincreases)
[taub tau h sig Z S sigma sen n senplot cilower ciupper] = ktaub([t;sigincreases]',0.3,1)
%We get a p-value of 0.087 if we average the three stations. Mayb we should
%fig out a way of doing this by weighting spatially so we could arrive at a
%"hotspot" of statistical certainty between these three stations, for
%example.

%54.53% decrease from the 76-85 baseline!

hold on
scatter(t,YearFreqrel(10,t-1975))
scatter(t,YearFreqrel(12,t-1975))
scatter(t,YearFreqrel(6,t-1975))
ylabel('Total Yearly Hrs. of FZRA')


%Finally for the weird increases on Long Island
% _________________________________________________________________________________________
%Try the mean or median station observation for each year instead:
t = 1976:2014
indices = find(a.StationLocations(:,1) > 40 & a.StationLocations(:,1) < 41.5 & ...
               a.StationLocations(:,2) < -72 & a.StationLocations(:,2) > -74) 
sigincreases = YearFreqrel(indices,t-1975)
sigincreases = mean(sigincreases)
[taub tau h sig Z S sigma sen n senplot cilower ciupper] = ktaub([t;sigincreases]',0.3,1)
%We get a p-value of 0.087 if we average the three stations. Mayb we should
%fig out a way of doing this by weighting spatially so we could arrive at a
%"hotspot" of statistical certainty between these three stations, for
%example.

%54.53% decrease from the 76-85 baseline!

hold on
scatter(t,YearFreqrel(10,t-1975))
scatter(t,YearFreqrel(12,t-1975))
scatter(t,YearFreqrel(6,t-1975))
ylabel('Total Yearly Hrs. of FZRA')


%% Boxplot   
%Years:
% startyear = 1976;
% endyear = 1996;
startyear = 2000;
endyear = 2014;
startyear = startyear - 1976;
endyear = endyear - 1976;

%Go through each station and create an average for each month for the time
%period chosen in "years":
for z = 1:97
    MonthFreqrel(z,1) = nanmean(a.MonthFreq_yearly_rel_series(z,((startyear*12)+1):12:(endyear*12)))/(31*24);
    MonthFreqrel(z,2) = nanmean(a.MonthFreq_yearly_rel_series(z,((startyear*12)+2):12:(endyear*12)))/((((28*29)+(29*10))/39)*24);
    MonthFreqrel(z,3) = nanmean(a.MonthFreq_yearly_rel_series(z,((startyear*12)+3):12:(endyear*12)))/(31*24);
    MonthFreqrel(z,4) = nanmean(a.MonthFreq_yearly_rel_series(z,((startyear*12)+4):12:(endyear*12)))/(30*24);
    MonthFreqrel(z,5) = nanmean(a.MonthFreq_yearly_rel_series(z,((startyear*12)+5):12:(endyear*12)))/(31*24);
    MonthFreqrel(z,6) = nanmean(a.MonthFreq_yearly_rel_series(z,((startyear*12)+6):12:(endyear*12)))/(30*24);
    MonthFreqrel(z,7) = nanmean(a.MonthFreq_yearly_rel_series(z,((startyear*12)+7):12:(endyear*12)))/(31*24);
    MonthFreqrel(z,8) = nanmean(a.MonthFreq_yearly_rel_series(z,((startyear*12)+8):12:(endyear*12)))/(31*24);
    MonthFreqrel(z,9) = nanmean(a.MonthFreq_yearly_rel_series(z,((startyear*12)+9):12:(endyear*12)))/(30*24);
    MonthFreqrel(z,10) = nanmean(a.MonthFreq_yearly_rel_series(z,((startyear*12)+10):12:(endyear*12)))/(31*24);
    MonthFreqrel(z,11) = nanmean(a.MonthFreq_yearly_rel_series(z,((startyear*12)+11):12:(endyear*12)))/(30*24);
    MonthFreqrel(z,12) = nanmean(a.MonthFreq_yearly_rel_series(z,((startyear*12)+12):12:(endyear*12)))/(31*24);
end


figure(19)
%boxplot([MonthFreqrel(:,10:12)*100 MonthFreqrel(:,1:5)*100])

%c = get(gca,'colororder') ; %save the color order for repetition of this plot
set(gca,'ColorOrder',c)
h = violinplot([MonthFreqrel(:,10:12)*100 MonthFreqrel(:,1:5)*100])
set(gca, 'XTick',1:8, 'XTickLabel',{'Oct' 'Nov' 'Dec' 'Jan' 'Feb' 'Mar' 'Apr' 'May'})
ylabel('Relative Frequency (%)')
grid on
hold on
line(1:8,median([MonthFreqrel(:,10:12)*100 MonthFreqrel(:,1:5)*100]))
set(gca,'FontSize',16)

%Save for comp.
MonthFreqrel1 = MonthFreqrel;

for month = 1:12
    [p(month) h(month)] = ranksum(MonthFreqrel1(:,month), MonthFreqrel(:,month))
end

%% Mann-Kendall Interpolation
t = 1976:2014;

for m = 1:97
    [taub(m) tau(m) h(m) sig(m) Z(m) S(m) sigma(m) sen(m) n senplot cilower(m) ciupper(m)] = ktaub([t;a.YearFreq_rel(m,t-1975)]',0.05,1);
end

%Percent change relative to 1976-1985 baseline:
 sen_percent = sen'./nanmean(YearFreqrel(:,1:10),2)*100*28;
 
%Set the boundaries to be near the local values:
sen_wboundpts = [sen'; mean(sen); (sen(75)+sen(73))/2; ...
    mean(sen); (sen(46)+sen(81))/2];
sig_wboundpts = [sig'; mean(sig); (sig(75)+sig(73))/2; ...
    mean(sig); (sig(46)+sig(81))/2];

senmap = griddata(stationlocations_wboundpts(:,1),stationlocations_wboundpts(:,2), ...
    sen_wboundpts,domain_x,domain_y,'natural'); %beautiful! Last parameter describes the method of interpolation.

sigmap = griddata(stationlocations_wboundpts(:,1),stationlocations_wboundpts(:,2), ...
    sig_wboundpts,domain_x,domain_y,'natural'); %beautiful! Last parameter describes the method of interpolation.
sigmap(find(sigmap > 0.1)) = NaN;

figure(11)
worldmap(lat_lim,lon_lim)
colormap(parula)
%caxis([-1, 1])
sen_ref = georasterref('RasterSize', size(senmap), 'Latlim', lat_lim, 'Lonlim', lon_lim);
geoshow(states,'FaceColor',[1 1 1])
hold on
geoshow(provinces,'FaceColor',[1 1 1])
framem('ffacecolor',[.5 .7 .9]); %shows water as blue
geoshow(senmap,sen_ref,'DisplayType','texturemap','FaceAlpha',0.9)
colorbar
%geoshow(sigmap,sen_ref,'DisplayType','texturemap','FaceAlpha',0.4)
plotm(a.StationLocations(:,1),a.StationLocations(:,2),'r+') %plots measurement stations

hold on
title('Linear Trend in Freezing Rain (2000-2014)')


%_________________________________________GOOD BUBBLE PLOT:
figure(1000)
worldmap(lat_lim,lon_lim)
caxis([-150,150])
cptcmap('SVS_tempanomaly', 'mapping', 'scaled');
c = colorbar; 

geoshow(states,'FaceColor',[.8 .8 .8])
hold on
geoshow(provinces,'FaceColor',[.8 .8 .8])
framem('ffacecolor',[.5 .7 .9]); %shows water as blue
%geoshow(fzramap,fzra_ref,'DisplayType','texturemap','FaceAlpha',0.8)
plotm(a.StationLocations(:,1),a.StationLocations(:,2),'r+') %plots measurement stations

sizecoeff = 700;
%scatterm(a.StationLocations(:,1),a.StationLocations(:,2),500*(1-sig),hoursperyear,'filled')
scatterm(a.StationLocations(:,1),a.StationLocations(:,2),500*(1-sig),sen_percent,'filled')
%Plot green Xs on locations that show significant trends
plotm(a.StationLocations(sig<0.05,1),a.StationLocations(sig<0.05,2),'gx') %plots measurement stations
c = colorbar
c.FontSize = 18
%c.FontName = 'Gill Sans MT'
ypercentages = get(c,'YTickLabel');
perc = repmat('%',size(ypercentages,1),1);
ypercentages = strcat(ypercentages, perc);
set(c,'YTickLabel',ypercentages);
xlabel(c,'Change from 1976-1985 Baseline')
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
                = sktt([tyears;tmonths;a.MonthFreq_yearly_rel_series(75,:)]',0.8,0,1)
            

%Alright, let's do this thing old-skool.
for m = 1:97
    tic
    for month = 1:12
        [taub(m,month) tau(m,month) h(m,month) sig(m,month) Z(m,month) S(m,month) sigma(m,month) sen(m,month) n senplot cilower(m,month) ciupper(m,month)] = ktaub([t;a.MonthFreq_yearly_rel_series(m,month:12:end)]',0.05,1);
    end
    toc
end

%Let's do it now after just taking the regional mean or median. And use
%this to apply to our seasonal frequencies:
for month = 1:12
    monthmean(month,:) = nanmean(a.MonthFreq_yearly_rel_series(:,month:12:end));
    %monthmean(month,:) = median(a.MonthFreq_yearly_rel_series(:,month:12:end),'omitnan');
    [taub(month) tau(month) h(month) sig(month) Z(month) S(month) sigma(month) sen(month) n senplot cilower(month) ciupper(month)] = ktaub([t;monthmean(month,:)]',0.05,1);
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
plotm(a.StationLocations(:,1),a.StationLocations(:,2),'r+') %plots measurement stations

hold on
title('Linear Trend in Freezing Rain (2000-2014)')
