%% Matt_FZRA_temps.m
% Updated October 2017

load b_fixedcoords.mat

%% Make a histogram of temperatures during FZRA.
stationnames = fieldnames(b);

temps = [];
coords = temps;
windspd = temps;
winddir = temps;
pctobs = temps;

for k = 1:length(stationnames)
    temps = [temps; b.(stationnames{k}).TEMP];
    coords = [coords; b.(stationnames{k}).coords];
    
    spd = b.(stationnames{k}).SPD;
    dir = b.(stationnames{k}).DIR;
    spd(spd == -999) = [];
    dir(dir == -999) = [];

    windspd = [windspd; mean(spd(1:ceil(0.36*length(spd))))];
    winddir = [winddir; mode(dir(1:ceil(0.36*length(spd))))];
    pctobs  = [pctobs;  mean(b.(stationnames{k}).numobsperYR)/8760];
end

%Eliminate -99 vals:
temps(temps == -999) = [];

figure(10)
hist(temps,-10:1:40)

stdev = std(convtemp(temps,'F','C'))
convtemp(mean(temps),'F','C')


%% Initial mapping and wind map
%Specify a map area of interest (the bounds of our map):
lat_lim = [38 50]; %deg N
lon_lim = [-95 -71]; %deg W

%Import the shapefile:
states = shaperead('usastatehi','UseGeoCoords',true);
provinces = shaperead('province','UseGeoCoords',true);

%Let's plot this along with the stations to make sure we've got the coastline situated correctly.
figure(2)
worldmap(lat_lim,lon_lim) %plots empty axes
title('Wind During FZRN Events in the Great Lakes');
geoshow(states,'FaceColor',[.5 1 .5])
geoshow(provinces,'FaceColor',[.5 1 .5])
framem('ffacecolor',[.5 .7 .9]); %shows water as blue
plotm(coords(:,1),coords(:,2),'r+') %plots measurement stations with red star
%Fix the wind dirs:
quiverm(coords(:,1),coords(:,2),cosd(-winddir).*windspd,sind(-winddir).*windspd,0.05);
xlabel('Longitude (deg W)')
ylabel('Latitude (deg N)')
title('')


%% Make a polar scatter plot of he wind directions for the four S Canada Stations:

figure(2)
worldmap(lat_lim,lon_lim) %plots empty axes
title('Wind During FZRN Events in the Great Lakes');
geoshow(states,'FaceColor',[.5 1 .5])
geoshow(provinces,'FaceColor',[.5 1 .5])
framem('ffacecolor',[.5 .7 .9]); %shows water as blue
plotm(coords(:,1),coords(:,2),'r+') %plots measurement stations with red star
quiverm(coords(:,1),coords(:,2),cosd(-winddir).*windspd,sind(-winddir).*windspd,0.05);

hold on

figure(3)
% figure(4)
for k = 1:97
    spd1 =[];
    spd = [];
    dir = [];
    spd1 = b.(stationnames{k}).SPD;
    dir = b.(stationnames{k}).DIR;
    spd = spd1;
    spd(spd1 == -999) = [];
    dir(spd1 == -999) = [];

    windspdmean = mean(spd(1:ceil(0.36*length(spd))));
    winddirmedian = median(dir(1:ceil(0.36*length(spd))));
    
    figure(3)
    polarhistogram(dir*180/pi,12)
    hold on
    polarplot([0,winddirmedian*180/pi],[0,50])
    
%     [X,Y] = pol2cart(dir*180/pi,spd);
%     figure(4)
%     hist3([X',Y'],[100,100])
    %polarscatter(dir*180/pi, spd)
    
    
    
    figure(2)
    plotm(coords(k,1),coords(k,2),'b*')
    
    children = get(gca, 'children');
    
    z = waitforbuttonpress;
    delete(children(10));
    
end
