% FZRA Intensity Investigation with CPC Precipitation Data
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

%wget --user=mairish --ask-#Gumpzilla -b -i ftp://ftp.cdc.noaa.gov/Datasets/cpc_us_precip/precip.V1.0.1976.nc
%To loop through files:
%for i in {77..99}; do wget --user=mairish --ask-#Gumpzilla -b -i ftp://ftp.cdc.noaa.gov/Datasets/cpc_us_precip/precip.V1.0.19${i}.nc; done

%wget -b -i ftp://ftp.cdc.noaa.gov/Datasets/cpc_us_precip/precip.V1.0.19${i}.nc

clear 

load a_all
load b
stationnames = fieldnames(b);

%Import station regional cluster identities, too:
load clusters_n_k5stationclusters_BETTER_plus_NYC_6thcluster

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

%Read in the CPC monthlies file:
%oid = 'cfrzr.1979.nc';
oid = 'precip.V1.0.mon.mean.nc'
ncdisp(oid);
nc_dump(oid) %SNCTOOLS version of ncdisp
setpref('SNCTOOLS','PRESERVE_FVD',true) %sets fastest-varying dimension first. Whether or not to transpose. It's a performance thing. - from http://mexcdf.sourceforge.net/tutorial/ch02.html
today = datestr(floor(now),'yyyymmdd'); %unused rmeow but this is a good reference command

f = ncinfo(oid);
nvars = length(f.Variables);
for k = 1:nvars
   varname= f.Variables(k).Name;
   disp(['Reading:  ' varname]);
   eval([varname ' = ncread(''' oid ''',''' varname ''');']);
   
end

%Convert from the "hours since 1700" to a datetime vector.
months = datetime(time*60*60,'ConvertFrom','epochtime','Epoch','1700-01-01');


%% Plot a monthly as a starter!! Take a look at January, 1976 for starters.
%Actually, we'll start by just averaging for the whole year.
lats = double(lat); %these are for the data we're plotting. Theyre singles by default from ESRL.
lons = double(lon);
precip = double(precip);
[domain_x, domain_y] = meshgrid(lons,lats);
latlim_precip = [min(lats) max(lats)];
lonlim_precip = [min(lons) max(lons)];

precipz = permute(precip,[2 1 3]);
precipz(precipz < 0) = NaN;

%Import the shapefile:
states = shaperead('usastatehi','UseGeoCoords',true);
provinces = shaperead('province','UseGeoCoords',true);

%Overlay the study region and CPC availability on one plot:
figure(2)
title('Jan 1976 Average')
worldmap(lat_lim,lon_lim) %use the one below this generally. not worldmap.
precip_ref = georasterref('RasterSize', size(precipz), 'Latlim', latlim_precip, 'Lonlim', lonlim_precip);
geoshow(states,'FaceColor',[1 1 1])
geoshow(provinces,'FaceColor',[1 1 1])
%framem('ffacecolor',[.5 .7 .9]); %shows water as blue
geoshow(precipz(:,:,1),precip_ref,'DisplayType','texturemap','FaceAlpha',0.8)
%plotm([sort([lat_lim_NARR lat_lim_NARR]) lat_lim_NARR(1)],[lon_lim_NARR fliplr(lon_lim_NARR) lon_lim_NARR(1)])
framem; gridm
caxis([0 40])
c = colorbar
c.FontSize = 22
c.FontName = 'Gill Sans MT'
tightmap %<---------------BEST FUNCTION EVER OMG


%% Parse the monthly precip values at each station and make a matrix of their trends over the years. 

%Now record the monthly precip value for each year for Nov. to April at
%each station, making a matrix that can be averaged for each region and
%plotted as time series!
monthstosample = months(months >= '01-Nov-1976 00:00:00' & months < '01-Jan-2015 00:00:00');
monthstosample(month(monthstosample) < 11 & month(monthstosample) > 4) = []; %get rid of May-Oct

%Find which grids are the closest to each station in the CPC dataset:
for ind = 1:length(stationnames)
    [c index] = min(abs(b.(stationnames{ind}).coords(1) - lats))
    closestCPCgridcoords(ind,1) = lats(index)
    closestCPCgridindices(ind,1) = index
    [c index] = min(abs(b.(stationnames{ind}).coords(2) - [lons-360]))
    closestCPCgridindices(ind,2) = index
end

%Great! Now let's actually build the matrix:
for sample = 1:length(monthstosample)
    for station = 1:length(stationnames)
        %mm/day average per month:
        precip_time_series(station,sample) = precipz(closestCPCgridindices(station,1),closestCPCgridindices(station,2),months == monthstosample(sample));
        if rem(sample,6)==0
            precip_time_series_yrly(station,sample/6) = mean(precip_time_series(station,sample-5:sample),'omitnan');
            precip_time_series_time(sample/6) = monthstosample(sample);
        end
    end
end

%And now let's average these time series into only six--one for each
%region!
for w = 1:6 %there are five clusters
    precip_time_series_regional_avg(w,:) = mean(precip_time_series(IDX_stations == w,:),'omitnan');
    precip_time_series_yrly_regional_avg(w,:) = mean(precip_time_series_yrly(IDX_stations == w,:),'omitnan');
    
end

%OMG the time hath come! Plot 'em and hope for the best!
figure(40)
for w = 1:6
    plot(precip_time_series_time,precip_time_series_yrly_regional_avg(w,:),'color',colorz(w,:),'LineWidth',2)
    hold on
    %Now calculate a least squares fit for the monthly time series and plot that, too:
    coeffs(w,:) = polyfit(datenum(monthstosample)',precip_time_series_regional_avg(w,:),1);
    bestfit(w,:) = polyval(coeffs(w,:),datenum(monthstosample)');
    corrc = corrcoef(bestfit(w,:),precip_time_series_regional_avg(w,:))
    rsquared(w) = corrc(1,2).^2;
    
    plot(monthstosample,bestfit(w,:),'color',colorz(w,:),'LineWidth',6)
        
    stats = fitlm(datenum(monthstosample)',precip_time_series_regional_avg(w,:),'linear');
    pvals(w) = table2array(stats.Coefficients(1,4));
end

xlabel('Year')
ylabel('Avg. Precipitation (mm/day)')
%legend('Northwest','Appalachia','Eastern Hotspots','S.-Central Canada','South Central','NYC')


%% Now let's do the same plot for overall trends in frequency for comparison:
%Let's average these time series into only six--one for each
%region!
t = 1976:2014;

for m = 1:97
    [taub(m) tau(m) h(m) sig(m) Z(m) S(m) sigma(m) sen(m) n senplot cilower(m) ciupper(m)] = ktaub([t;a.YearFreq_rel(m,t-1975)]',0.05,1);
end


for w = 1:6 %there are five clusters
    freq_regional_avg(w,:) = mean(a.YearFreq(IDX_stations == w,:),'omitnan');
    freq_sen_regional_avg(w) = mean(sen(IDX_stations == w),'omitnan');
    
end

%OMG the time hath come! Plot 'em and hope for the best!
figure(41)
for w = 1:6
    plot(1976:2014,freq_regional_avg(w,:),'color',colorz(w,:),'LineWidth',2)
    hold on
    %Now calculate a least squares fit for the monthly time series and plot that, too:
    coeffs_freq(w,:) = polyfit(1976:2014,freq_regional_avg(w,:),1);
    bestfit_freq(w,:) = polyval(coeffs_freq(w,:),1976:2014);
    corrc = corrcoef(bestfit_freq(w,:),freq_regional_avg(w,:))
    rsquared_freq(w) = corrc(1,2).^2;
    
    plot(1976:2014,bestfit_freq(w,:),'color',colorz(w,:),'LineWidth',6)
    
    stats = fitlm(1976:2014,freq_regional_avg(w,:),'linear');
    pvals_freq(w) = table2array(stats.Coefficients(1,4));
    
end

xlabel('Year')
ylabel('Hours of FZRA per Year')
%legend('Northwest','Appalachia','Eastern Hotspots','S.-Central Canada','South Central','NYC')



    
    


