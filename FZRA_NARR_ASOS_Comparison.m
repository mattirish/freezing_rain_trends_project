% FZRA NARR-ASOS Comparison
%
% Expand the path to the entire MATLAB folder before running this to access
% SNCTOOLS, the package I downloaded that ESRL recommends for working with
% its data.
% Great SNCTOOLS tutorial: http://mexcdf.sourceforge.net/tutorial/
% NARR model data is output onto a 349x277 Lambert Conformal Conic grid
% (called "NCEP Grid 221")
%
% See ESRL PSD's catalog here for OPeNDAP file addresses!
% https://www.esrl.noaa.gov/psd/thredds/catalog/Datasets/NARR/catalog.html
%
% WHENEVER YOUR MAP LAYOUTS BRING YOU WOES, JUST USE TIGHTMAP!!! <3 tightmap <3
%
% Matt Irish 2018

clear 

% set(0,'DefaultAxesFontSize', 15 );
% set(0,'DefaultAxesFontWeight', 'bold');
% set(0,'DefaultTextFontSize',15);
% set(0,'DefaultTextFontWeight','bold');
% %set(0,'DefaultLineLineWidth',4);

%% Open dis puppy.
%Boundaries of our study region:
lat_lim = [38 50]; %deg N
lon_lim = [-95 -71]; %deg W
lat_lim_narrdef = [0 80]; %deg N default NARR plot
lon_lim_narrdef = [150 0]; %deg W
lat_lim_NARR = [25 55];
lon_lim_NARR = [-100 -65];

ncdisp('cfrzr.1979.nc')

%Read in files from the interwebs through OPeNDAP.
%URL below is daily averages: 
url = 'http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/NARR/Dailies/pressure/air.197901.nc';
%Full 3-hrly file for January 1979:
url = 'http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/NARR/pressure/air.197901.nc';

nc_dump(url)

nc_dump('cfrzr.1979.nc') %SNCTOOLS version of ncdisp
setpref('SNCTOOLS','PRESERVE_FVD',true) %sets fastest-varying dimension first. Whether or not to transpose. It's a performance thing. - from http://mexcdf.sourceforge.net/tutorial/ch02.html
today = datestr(floor(now),'yyyymmdd'); %unused rmeow but this is a good reference command

oid='cfrzr.1979.nc';

f = ncinfo(oid);
nvars = length(f.Variables);
for k = 1:nvars
   varname= f.Variables(k).Name;
   disp(['Reading:  ' varname]);
   eval([varname ' = ncread(''' oid ''',''' varname ''');']);
end


%Now find limits in terms of x and y, the variables from the NARR, which
%are in meters :( one degree lat is about 111 km, one deg lon is about 100
%km
% [dum yinds] = min(abs(2775000-y)) %southern lat bound
% [dum yindn] = min(abs(6105000-y)) %northern lat bound
% yindn = yindn - yinds;  %counting offset
% [dum xindw] = min(abs(5565600-x)) %western lon bound
% [dum xinde] = min(abs(8580750-x)) %eastern lon bound
% xinde = xinde - xindw;  %counting offset
[dum yinds] = min(abs(2305000-y)); %southern lat bound
[dum yindn] = min(abs(5705000-y)); %northern lat bound
yindn = yindn - yinds;  %counting offset
[dum xindw] = min(abs(5865600-x)); %western lon bound
[dum xinde] = min(abs(9080750-x)); %eastern lon bound
xinde = xinde - xindw;  %counting offset

%Just use nc_varget to get a subset of the domain! We need this to make
%sure our lat and lon matrices are correct.

time = nc_varget('cfrzr.1979.nc','time'); %time
lat = nc_varget('cfrzr.1979.nc','lat',[yinds xindw],[yindn xinde]);
lon = nc_varget('cfrzr.1979.nc','lon',[yinds xindw],[yindn xinde]);
y = nc_varget('cfrzr.1979.nc','y',[yinds],[yindn]);
x = nc_varget('cfrzr.1979.nc','x',[xindw],[xinde]);
cfrzr = nc_varget('cfrzr.1979.nc','cfrzr',[0 yinds xindw],[inf yindn xinde]); %x,y,time


%% Plot that pupper! Take a look at Feb 23, 1979 for starters. Should see something.
%Actually, we'll start by just averaging for the whole year.
lats = double(lat); %these are for the data we're plotting. Theyre singles by default from ESRL.
lons = double(lon);

%Import the shapefile:
states = shaperead('usastatehi','UseGeoCoords',true);
provinces = shaperead('province','UseGeoCoords',true);

%Make a longitude matrix that flips deg East records to negative for
%plotting without the Mapping Toolbox:
lons4regplots = lons;
lons4regplots(lons > 100) = lons(lons>100)-360;

%Try and plot without Mapping Toolbox:
figure(1)
pcolor(lons4regplots,lats,squeeze(cfrzr(1,:,:)));shading interp;
load coast;hold on;plot(long,lat);plot(long+360,lat);hold off
colorbar


%Plot yearly total CFRZR using the Mapping Toolbox:
fzra_mean = squeeze(sum(cfrzr,1));
fzra_mean(fzra_mean < 0) = NaN;

%Overlay the study region, PCA region (bigger) and NARR availability on one plot:
figure(2)
worldmap([15 80],[-140 -30]) %use the one below this generally. not worldmap.
% axesm('eqdcylin','maplatlimit',lat_lim,'maplonlimit',lon_lim,...
%     'PLineLocation',5,'MLineLocation',5,...
%     'ParallelLabel','on','MeridianLabel','on','MLabelParallel','south') %plots empty axes

%axesm('eqdcylin','maplatlimit',[-80 80],'maplonlimit',[0 360]);  % Create a cylindrical equidistant map
geoshow(states,'FaceColor',[1 1 1])
geoshow(provinces,'FaceColor',[1 1 1])
%framem('ffacecolor',[.5 .7 .9]); %shows water as blue
pcolorm(lats,lons,fzra_mean,'FaceAlpha',0.8)           % pseudocolor plot "stretched" to the grid
plotm([sort([lat_lim lat_lim]) lat_lim(1)],[lon_lim fliplr(lon_lim) lon_lim(1)])
%plotm([sort([lat_lim_NARR lat_lim_NARR]) lat_lim_NARR(1)],[lon_lim_NARR fliplr(lon_lim_NARR) lon_lim_NARR(1)])
framem; gridm
caxis([0 40])
c = colorbar
c.FontSize = 22
c.FontName = 'Gill Sans MT'
tightmap %<---------------BEST FUNCTION EVER OMG


%% Read in data using OPeNDAP to recreate Cortinas 2000, Fig. 10. 
%First we'll just do a big old 2013 ice storm that formed a band from
%west of Chi across the Hudson River Valley & hit Toronto hard.
%Looks like it occurred from 21 Dec ~19:00 to 22 Dec ~8:00 2013.
daytoplot = '21-Dec-2013 21:00:00'; %only works for 8-hrly times rmeow, would need to add script to round to nearest 3-hrly record


%Create the url for each variable we want. Pressure for most and monolevel for cfrzr, etc:
%Air temps (K):
url_air = ['http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/NARR/pressure/air.',datestr(daytoplot,'yyyy'),datestr(daytoplot,'mm'),'.nc'];
%Geopotential height (m):
url_hgt = ['http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/NARR/pressure/hgt.',datestr(daytoplot,'yyyy'),datestr(daytoplot,'mm'),'.nc'];
%Winds (m/s):
url_uwnd = ['http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/NARR/pressure/uwnd.',datestr(daytoplot,'yyyy'),datestr(daytoplot,'mm'),'.nc'];
url_vwnd = ['http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/NARR/pressure/vwnd.',datestr(daytoplot,'yyyy'),datestr(daytoplot,'mm'),'.nc'];

%Monolevel stuff:
%Mean sea level pressure (Pa):
url_prmsl = ['http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/NARR/monolevel/prmsl.',datestr(daytoplot,'yyyy'),'.nc'];
%Relative Humidity at 2m (%):
url_rhum = ['http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/NARR/monolevel/rhum.2m.',datestr(daytoplot,'yyyy'),'.nc'];


%Figure out what the time index we want is:
timez = nc_varget(url_air,'time');
timez = datetime(timez*3600,'ConvertFrom','epochtime','Epoch','1800-01-01');
indextoplot = find(datetime(daytoplot) == timez) - 1; %minus 1 bc nc_varget uses zero-based indexing

%Download that data! Temps, heights, and winds at relevant levels:
tic; tic
air.at850 = nc_varget(url_air,'air',[indextoplot 6 xindw yinds],[1 1 xinde-xindw yindn-yinds]) - 273.15; %format is (x,y,level,time)
%air.at850 = nc_varget(url_air,'air',[0 0 6 indextoplot],[inf inf 1 1]) - 273.15; %format is (x,y,level,time)
toc

%%%%%%%The order of variables in all the below need to be fixed. it's time,level,x,y, not x,y,level,time. 
% air.at700 = nc_varget(url_air,'air',[20 29 12 indextoplot],[1 1 1 1]) - 273.15; %format is (time,level,x,y)
% air.at500 = nc_varget(url_air,'air',[0 0 16 indextoplot],[inf inf 1 1]) - 273.15;
% 
% hgt.at850 = nc_varget(url_hgt,'hgt',[0 0 6 indextoplot],[inf inf 1 1]);
% hgt.at700 = nc_varget(url_hgt,'hgt',[0 0 12 indextoplot],[inf inf 1 1]); 
% hgt.at500 = nc_varget(url_hgt,'hgt',[0 0 16 indextoplot],[inf inf 1 1]); 
% 
% uwnd.at850 = nc_varget(url_uwnd,'uwnd',[0 0 6 indextoplot],[inf inf 1 1]);
% uwnd.at700 = nc_varget(url_uwnd,'uwnd',[0 0 12 indextoplot],[inf inf 1 1]);
% uwnd.at500 = nc_varget(url_uwnd,'uwnd',[0 0 16 indextoplot],[inf inf 1 1]); 
% 
% vwnd.at850 = nc_varget(url_vwnd,'vwnd',[0 0 6 indextoplot],[inf inf 1 1]); 
% vwnd.at700 = nc_varget(url_vwnd,'vwnd',[0 0 12 indextoplot],[inf inf 1 1]);
% vwnd.at500 = nc_varget(url_vwnd,'vwnd',[0 0 16 indextoplot],[inf inf 1 1]);
toc
%Read in the monolevel stuff:
timez = nc_varget(url_prmsl,'time'); %Figure out new index.
timez = datetime(timez*3600,'ConvertFrom','epochtime','Epoch','1800-01-01');
indextoplot = find(datetime(daytoplot) == timez) - 1;

% prmsl = nc_varget(url_prmsl,'prmsl',[0 0 indextoplot],[inf inf 1]); %format is (x,y,time)
% rhum = nc_varget(url_rhum,'rhum',[0 0 indextoplot],[inf inf 1]); %format is (x,y,time)
prmsl = nc_varget(url_prmsl,'prmsl',[0 0 5],[inf inf 1]); %format is (x,y,time)
rhum = nc_varget(url_rhum,'rhum',[0 0 5],[inf inf 1]); %format is (x,y,time)


%Now read in monthly averages so we can derive anomaly maps from them!
%Figure out what the time index we want is:
timez_monthly = nc_varget(url_prmsl_monthly,'time');
timez_monthly = datetime(timez_monthly*3600,'ConvertFrom','epochtime','Epoch','1800-01-01');
%Data go from Jan. 1979 to Dec. 2016, and we only want until Dec. 2014 (up to and including index 432).

%Mean sea level pressure (Pa):
load monthlies.mat      %load saved giant monthly 3x3 arrays only when manipulating them here.
url_prmsl_monthly = ['http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/NARR/Monthlies/monolevel/prmsl.mon.mean.nc'];
%prmsl_monthly = nc_varget(url_prmsl_monthly,'prmsl',[0 0 0],[inf inf 432]); %format is (x,y,time)
prmsl_monthly_trim = prmsl_monthly(:,:,1:432); 
%Add 1000-500 hPa geopotential heights later, in accordance with Erfani et al.
%2012.

%Now make a 3x3 monthly averages array where [lat lon month].
%Each map is the 1979-2014 average for that month.
for z = 1:12
    prmsl_monthly_avg(:,:,z) = mean(prmsl_monthly_trim(:,:,z:12:end),3);
end

%% Do the same now for the 500 mb avg monthly heights:
url_hgtavg = 'http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/NARR/Monthlies/pressure/hgt.mon.mean.nc';
hgt500_monthly_avg = nc_varget(url_hgtavg,'hgt',[0 16 yinds xindw],[432 1 yindn xinde]);
hgt750_monthly_avg = nc_varget(url_hgtavg,'hgt',[0 10 yinds xindw],[432 1 yindn xinde]);
hgt850_monthly_avg = nc_varget(url_hgtavg,'hgt',[0 6 yinds xindw],[432 1 yindn xinde]);
hgt1000_monthly_avg = nc_varget(url_hgtavg,'hgt',[0 0 yinds xindw],[432 1 yindn xinde]);

hgt500_monthly = squeeze(hgt500_monthly_avg);
hgt750_monthly = squeeze(hgt750_monthly_avg);
hgt850_monthly = squeeze(hgt850_monthly_avg);
hgt1000_monthly = squeeze(hgt1000_monthly_avg);

%Now clear the monthly_avgs before running this loop. It needs to overwrite
%them.
%Now make an average of map.
for z = 1:12
    hgt500_monthly_avg(z,:,:) = mean(hgt500_monthly(z:12:end,:,:),1);
    hgt750_monthly_avg(z,:,:) = mean(hgt750_monthly(z:12:end,:,:),1);
    hgt850_monthly_avg(z,:,:) = mean(hgt850_monthly(z:12:end,:,:),1);
    hgt1000_monthly_avg(z,:,:) = mean(hgt1000_monthly(z:12:end,:,:),1);
    z
end
%Make sure the dimensions are in the right order before using this in
%SynopticWeatherTyping.m. We have time,y,x and we need y,x,time:
hgt500_monthly_avg = permute(hgt500_monthly_avg,[2,3,1]);
hgt750_monthly_avg = permute(hgt750_monthly_avg,[2,3,1]);
hgt850_monthly_avg = permute(hgt850_monthly_avg,[2,3,1]);
hgt1000_monthly_avg = permute(hgt1000_monthly_avg,[2,3,1]);

%Save 'em!

%Plot these monthly average maps to ensure everything's right:
%Map boundaries. Make 'em wide
lat_lim = [35 55]; %deg N
lon_lim = [-98 -68]; %deg W

figure(300)
axesm('eqdcylin','maplatlimit',lat_lim,'maplonlimit',lon_lim,...
    'PLineLocation',5,'MLineLocation',5,...
    'ParallelLabel','on','MeridianLabel','on','MLabelParallel','south'); %plots empty axes
framem; gridm;
geoshow(states,'FaceColor',[1 1 1]);
geoshow(provinces,'FaceColor',[1 1 1]);
%framem('ffacecolor',[.5 .7 .9]); %shows water as blue
%pcolorm(double(lats),double(lons),fzra_mean,'FaceAlpha',0.8)           % pseudocolor plot "stretched" to the grid

prmsl_monthly_avg(prmsl_monthly_avg <0) = NaN; prmsl_monthly_avg_mb = prmsl_monthly_avg/100;
%pcolorm(double(lats),double(lons),prmsl,'FaceAlpha',0.8)
rhum(rhum <0) = NaN;
prmsl_labels = [900:5:1100]/10; %kPa

%Plot MSL pressure and RH:
[C h] = contourm(lats,lons,prmsl_monthly_avg_mb(:,:,12)/10,prmsl_labels,'LineWidth',2);
tightmap; 
colorbar















%% Now make Cortinas's plot!
%Map boundaries. Make 'em wide
lat_lim = [35 55]; %deg N
lon_lim = [-98 -68]; %deg W

figure(3)
subplot(2,2,1)
axesm('eqdcylin','maplatlimit',lat_lim,'maplonlimit',lon_lim,...
    'PLineLocation',5,'MLineLocation',5,...
    'ParallelLabel','on','MeridianLabel','on','MLabelParallel','south'); %plots empty axes
framem; gridm;
geoshow(states,'FaceColor',[1 1 1]);
geoshow(provinces,'FaceColor',[1 1 1]);
%framem('ffacecolor',[.5 .7 .9]); %shows water as blue
%pcolorm(double(lats),double(lons),fzra_mean,'FaceAlpha',0.8)           % pseudocolor plot "stretched" to the grid

prmsl(prmsl <0) = NaN; prmsl_mb = prmsl/100;
%pcolorm(double(lats),double(lons),prmsl,'FaceAlpha',0.8)
rhum(rhum <0) = NaN;
prmsl_labels = [900:5:1100]/10; %kPa
rhum_labels = [0:10:100]; % percent

%Plot MSL pressure and RH:
[Ca ha] = contourm(lats,lons,rhum,rhum_labels,'LineWidth',2,'LineColor','k','LineStyle',':');
[C h] = contourm(lats,lons,prmsl_mb/10,prmsl_labels,'LineWidth',2);
tightmap; 

%Figure a better way out to label RH. Too time intensive and close
%together. Maybe need to smooth the map for non-averaged maps.
%clabelm(Ca);

%Add an "L" over the low pressure system:
prmsl_window = prmsl_mb(lats > lat_lim(1) & lats < lat_lim(2) ...
    & lons > lon_lim(1) & lons < lon_lim(2)); %selects only inside our region
windowmin = min(prmsl_window);
lattext = lats(prmsl_mb == windowmin);
lontext = lons(prmsl_mb == windowmin);
textindex = find(lattext > lat_lim(1) & lattext < lat_lim(2) ...
    & lontext > lon_lim(1) & lontext < lon_lim(2));
textm(lattext(textindex),lontext(textindex), 'L', 'FontWeight','bold','FontSize',20,'Color','r')

%Add an "H" over any high pressure system (cutoff at 1020 mb)
windowmax = max(prmsl_window);
if windowmax > 1020 %mb
    lattext = lats(prmsl_mb == max(prmsl_window));
    lontext = lons(prmsl_mb == max(prmsl_window));
    textindex = find(lattext > lat_lim(1) & lattext < lat_lim(2) ...
        & lontext > lon_lim(1) & lontext < lon_lim(2));
    textm(lattext(textindex),lontext(textindex), 'H', 'FontWeight','bold','FontSize',20,'Color','b')
end

caxis([windowmin windowmax]/10) %kPa. Sets color ramp to the range of our region.
c = colorbar;
xlabel(c,'MSL Pressure (kPa)')

%------------------------------------------------------------
subplot(2,2,2)
axesm('eqdcylin','maplatlimit',lat_lim,'maplonlimit',lon_lim,...
    'PLineLocation',5,'MLineLocation',5,...
    'ParallelLabel','on','MeridianLabel','on','MLabelParallel','south') %plots empty axes
%worldmap(lat_lim,lon_lim)
framem; gridm
geoshow(states,'FaceColor',[1 1 1])
geoshow(provinces,'FaceColor',[1 1 1])
%framem('ffacecolor',[.5 .7 .9]); %shows water as blue
%pcolorm(double(lats),double(lons),fzra_mean,'FaceAlpha',0.8)           % pseudocolor plot "stretched" to the grid
air.at850(air.at850 <-100) = NaN;
hgt.at850(hgt.at850 < 0) = NaN;
% [c,h] = contourm(lats,lons,air.at850) 
% ch = get(h,'child'); alpha(ch,0.8)
% hold on
% contourm(lats,lons,hgt.at850,'k--')

hgtat850_labels = [1000:50:2000]/1000; %km
airat850_labels = [-20:5:15]; %deg C

%Plot 850mb temps and geopotential height:
[Cc hc] = contourm(lats,lons,hgt.at850/1000,hgtat850_labels,'LineWidth',2,'LineColor','k','LineStyle',':');
[Cd hd] = contourfm(lats,lons,air.at850,airat850_labels);
tightmap; 

%Figure a better way out to label RH. Too time intensive and close
%together. Maybe need to smooth the map for non-averaged maps.
clabelm(Cc);

%Add an "L" over the low pressure system:
hgtat850_window = hgt.at850(lats > lat_lim(1) & lats < lat_lim(2) ...
    & lons > lon_lim(1) & lons < lon_lim(2)); %selects only inside our region
windowmin = min(hgtat850_window);
lattext = lats(hgt.at850 == windowmin);
lontext = lons(hgt.at850 == windowmin);
textindex = find(lattext > lat_lim(1) & lattext < lat_lim(2) ...
    & lontext > lon_lim(1) & lontext < lon_lim(2));
textindex = textindex(1); %just take one element
textm(lattext(textindex),lontext(textindex), 'L', 'FontWeight','bold','FontSize',20,'Color','r')

%Add an "H" over any high pressure system (cutoff at 1020 mb)
windowmax = max(hgtat850_window);
% if windowmax > 1020 %mb
%     lattext = lats(prmsl_mb == max(prmsl_window));
%     lontext = lons(prmsl_mb == max(prmsl_window));
%     textindex = find(lattext > lat_lim(1) & lattext < lat_lim(2) ...
%         & lontext > lon_lim(1) & lontext < lon_lim(2));
%     textm(lattext(textindex),lontext(textindex), 'H', 'FontWeight','bold','FontSize',20,'Color','b')
% end

%****Fix to get windowmin and max for temp, not height.
%caxis([windowmin windowmax]) %K. Sets color ramp to the range of our region.
c = colorbar;
xlabel(c,'Temperature (deg C)')


%-------------------------------------------------------------
subplot(2,2,3)
axesm('eqdcylin','maplatlimit',lat_lim,'maplonlimit',lon_lim,...
    'PLineLocation',5,'MLineLocation',5,...
    'ParallelLabel','on','MeridianLabel','on','MLabelParallel','south') %plots empty axes
framem; gridm
geoshow(states,'FaceColor',[1 1 1])
geoshow(provinces,'FaceColor',[1 1 1])
%framem('ffacecolor',[.5 .7 .9]); %shows water as blue
%pcolorm(double(lats),double(lons),fzra_mean,'FaceAlpha',0.8)           % pseudocolor plot "stretched" to the grid
air.at700(air.at700 <-100) = NaN;
hgt.at700(hgt.at700 < 0) = NaN;
[c,h] = contourfm(double(lats),double(lons),air.at700) 
ch = get(h,'child'); alpha(ch,0.8)
hold on
contourm(double(lats),double(lons),hgt.at700,'k--')
tightmap

%------------------------------------------------------------
subplot(2,2,4)
% axesm('eqdcylin','maplatlimit',lat_lim,'maplonlimit',lon_lim,...
%     'PLineLocation',5,'MLineLocation',5,...
%     'ParallelLabel','on','MeridianLabel','on','MLabelParallel','south') %plots empty axes
worldmap(lat_lim,lon_lim)
%framem; gridm
geoshow(states,'FaceColor',[1 1 1])
geoshow(provinces,'FaceColor',[1 1 1])
%framem('ffacecolor',[.5 .7 .9]); %shows water as blue
%pcolorm(double(lats),double(lons),fzra_mean,'FaceAlpha',0.8)           % pseudocolor plot "stretched" to the grid
hgt.at500(hgt.at500 < 0) = NaN;
[handle h ] = contourm(double(lats),double(lons),hgt.at850,'k--')
clabelm(handle)
uwnd500plot = NaN.*zeros(size(uwnd.at500)); %we want to only plot one vector for every few cells
uwnd500plot(1:6:end,1:6:end) = uwnd.at500(1:6:end,1:6:end);
vwnd500plot = NaN.*zeros(size(vwnd.at500));
vwnd500plot(1:6:end,1:6:end) = vwnd.at500(1:6:end,1:6:end);
quiverm(double(lats),double(lons),vwnd500plot,uwnd500plot,10);
%quiverm(double(lats),double(lons),vwnd.at500,uwnd.at500,10);
colorbar

%Add a supertitle (has to be done last):
set(gcf,'NextPlot','add');
axes;
h = title(['Meteorological Conditions on ',daytoplot]);
set(gca,'Visible','off');
set(h,'Visible','on');



figure(1000)
worldmap() %plots empty axes
pcolorm(double(lats),double(lons),sqrt((vwnd.at500).^2+(uwnd.at500).^2),'FaceAlpha',.7)
caxis([0 40])

%should get this to make windbarbs https://www.mathworks.com/matlabcentral/fileexchange/33851-wind-barb-plotter


