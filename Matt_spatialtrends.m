%% Matt_spatialanalysis.m
% February 2016

load spatialtrends.mat
load frachoursavail.mat
load Wind.mat
load totalobsperstationperyear.mat

%Specify a map area of interest (the bounds of our map):
lat_lim = [41 48]; %deg N
lon_lim = [-92 -82]; %deg W
%Import the shapefile:
% coastline = shaperead('Lake_Michigan_Coast.shp','UseGeoCoords', true);
states = shaperead('usastatehi','UseGeoCoords',true);
provinces = shaperead('province','UseGeoCoords',true);

%Let's plot this along with the stations to make sure we've got the coastline situated correctly.
figure(1)
worldmap(lat_lim,lon_lim) %plots empty axes
title('Wind Measurement Stations on Lake Michigan');
geoshow(states,'FaceColor',[.5 1 .5])
geoshow(provinces,'FaceColor',[.5 1 .5])
framem('ffacecolor',[.5 .7 .9]); %shows water as blue
plotm(stationlocations(:,1),stationlocations(:,2),'r+') %plots measurement stations with red star
quiverm(stationlocations(:,1)',stationlocations(:,2)',cosd(-WindDir).*WindSpeed,sind(-WindDir).*WindSpeed,0.3);
xlabel('Longitude (deg W)')
ylabel('Latitude (deg N)')
northarrow('latitude', 41.5, 'longitude', -89.5) %adding north arrow for pizzazz

%% 2D-Interpolating between points to make a raster map of average freezing rain freq
%We can use the griddata function to interpolate the data using
%the cubic method (better than linear interpolation for modeling
%things that vary smoothly across space)
%Initialize domain variables, etc.
resolution = 0.01;
x = (lon_lim(1):resolution:lon_lim(2));
y = (lat_lim(1):resolution:lat_lim(2));
[domain_y, domain_x] = meshgrid(x,y);

%Throw out Benton Harbor and Iron Mountain:
% stationlocations(8,:) = [];
% stationlocations(19,:) = [];
% YearFreq(8,:) = [];
% YearFreq(19,:) = [];
% MonthFreq(8,:) = [];
% MonthFreq(19,:) = [];
% frachoursavail(8) = [];
% frachoursavail(19) = [];

%Normalize to the total number of hours reported at each station:
numobservations = frachoursavail*350264; %scale fracs by total hours btwn 1976 and end of 2014
YearFreqrel = bsxfun(@times,YearFreq,1./(numobservations./39)); %divide each col of YearFreq by vec of est total num observations
YearFreqrel = YearFreqrel*100*100; %scale by 100, then make percent

%Also include the four corner points of our region of interest with average
%frequencies so that the interpolation is carried out on the whole region:
stationlocations_wboundpts = [stationlocations; lat_lim' lon_lim'; lat_lim' fliplr(lon_lim)'];
MonthFreq_wboundpts = [MonthFreq; ones(1,12).*mean(MonthFreq); ones(1,12).*mean(MonthFreq); ...
    ones(1,12).*mean(MonthFreq); ones(1,12).*mean(MonthFreq)];
YearFreqrel_wboundpts = [YearFreqrel; ones(1,39).*mean(YearFreqrel); ones(1,39).*mean(YearFreqrel); ...
    ones(1,39).*mean(YearFreqrel); ones(1,39).*mean(YearFreqrel)];

%It's time! Interpolate the scattered data points using the griddata func:
yearstoplot = 1990:2014; %select range of years to map average of
fzrnvals_wboundpts = mean(YearFreqrel_wboundpts(:,yearstoplot - 1975),2);
fzrnmap = griddata(stationlocations_wboundpts(:,1),stationlocations_wboundpts(:,2),fzrnvals_wboundpts,domain_x,domain_y,'natural'); %beautiful! Last parameter describes the method of interpolation.


% fzrnmap19902004 = griddata(stationlocations_wboundpts(:,1),stationlocations_wboundpts(:,2),fzrnvals_wboundpts,domain_x,domain_y,'natural'); %beautiful! Last parameter describes the method of interpolation.
% fzrnmapdiff = fzrnmap - fzrnmap19902004;

% %Now we need a boundary matrix to mask the map with.
% %Convert the vector boundaries of US to a matrix using the vec2mtx function:
% cellsperdegree = 1/resolution; %resolution of map raster
% [boundary_matrix, R_dummy] = vec2mtx(states(22).Lat, states(22).Lon, cellsperdegree, lat_lim, lon_lim, 'filled');
% 
% %Add one row and column to the US boundary matrix with its "out of
% %bounds" default value of 2 to make it match the size of the raster matrix
% us_boundary = 2*ones(size(fzrnmap));
% us_boundary(2:end,2:end) = boundary_matrix;

%Plot the masked windmap! She's gonna be beautiful.
figure(2)
worldmap(lat_lim,lon_lim)
colormap(parula)
demcmap(fzrnmap,64,'window','window') %applies colormap to the raster data
fzrn_ref = georasterref('RasterSize', size(fzrnmap), 'Latlim', lat_lim, 'Lonlim', lon_lim);
% geoshow(states,'FaceColor',[.5 1 .5])
% geoshow(provinces,'FaceColor',[.5 1 .5])
geoshow(states,'FaceColor',[1 1 1])
geoshow(provinces,'FaceColor',[1 1 1])
framem('ffacecolor',[.5 .7 .9]); %shows water as blue
geoshow(fzrnmap,fzrn_ref,'DisplayType','texturemap','FaceAlpha',0.8)
plotm(stationlocations(:,1),stationlocations(:,2),'r+') %plots measurement stations
colorbar
% geoshow(fzrnmap,fzrn_ref,'DisplayType','contour')
[C h] = contourm(fzrnmap,fzrn_ref,'LevelStep',5,'LineWidth',1,'LineColor','k');
clabelm(C)
title('Mean Yearly Frequency of Hourly FZRN Reports (1976 - 1990)')
northarrow('latitude', 41.5, 'longitude', -89.5)

%% boxplot
for k = 1: length(MonthFreq)
MonthFreq_rel(k,:) = MonthFreq(k,:)./mean(totalobsperstationperyear(k,:));
end
MonthFreq_rel(8,:) = [];
MonthFreq_rel(11,:) = [];
MonthFreq_rel(17,:) = [];

figure(3)
boxplot(MonthFreq_rel*100)
xlabel('Month')
ylabel('Frequency (%)')




