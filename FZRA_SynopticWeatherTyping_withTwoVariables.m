% FZRA Automated Synoptic Weather Typing
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
% WHENEVER YOUR MAP LAYOUTS BRING YOU WOES, JUST USE TIGHTMAP!!! 
%
% Matt Irish 2018

clear 

%load monthlies_correctedfordomainsize %archive of MSL maps for all 456 months and 12 averages, one for each month!

%Specify a map area of interest (the bounds of our map):
lat_lim = [37.5 50]; %deg N
lon_lim = [-95 -71]; %deg W

%% Load the system boundaries.
%Load the lat, lon, x, and y, (along with cfrzr, from 1979), the subregion
%domain I chose. Routine is at the beginning of
%FZRA_SynopticWeatherTyping.m.
load latslons_subregion_for_PCA


%% Reduce the domain of the monthlies down to our chosen region.
%prmsl_monthly_avg = prmsl_monthly_avg(xindw+1:xindw+xinde,yinds+1:yinds+yindn,:);
% %Rearrange so that it goes y,x,time:
%prmsl_monthly_avg = permute(prmsl_monthly_avg,[2,1,3]);
%Okay, I think we're good to go now. The dim order is y by x by time.

%% Download surface MSL maps for all the events of 2014 and create anomaly 
% maps of each by subtracting the event average from the all-time monthly
% average.

%Call FZRA_EventTimes to give us a 3-hourly NARR-ready log of all events:
timestep = 3;               %rounds to every three hours
min_reports_per_event = 1;  %min no. reports that constitute an event. 4 is the saved files.
max_nonevent_hrs = 6;       %allow up to 4 hours between events

%[event_times event_ids] = FZRA_EventTimes(timestep, min_reports_per_event, max_nonevent_hrs);
[event_times, event_ids, event_spd, event_spd_std, event_precip, event_precip_std, event_pct_lightFZRA, event_stationcounts, nonevent_times, nonevent_stations] = FZRA_EventTimes(timestep, min_reports_per_event, max_nonevent_hrs);
%load allevents_ids_&_times (This is 343. Idk if it's any different than eventtimesoutput_343. I think it's not.)
%load eventtimesoutput_343 %(343)

%load eventtimesoutput_311 %(311) basically, any report is an event and there's no grace period between reports.
%load eventtimesoutput_316 %(316) basically, any report is an event and there's a big grace period between reports.

%Select only dates for the decided time window:
% event_ids = event_ids(event_times >= datetime(1997,1,1) & event_times < datetime(2015,1,1));

event_times_case = event_times(event_times >= datetime(1979,1,1) & event_times < datetime(2015,1,1));
%event_times_case = event_times(event_times >= datetime(1997,1,1) & event_times < datetime(2015,1,1));
%event_ids = event_ids(event_times >= datetime(1979,1,1) & event_times < datetime(1997,1,1));

%event_times_case = event_times(event_times >= datetime(1979,1,1) & event_times < datetime(1997,1,1));

%Run the following loop if you're not downloading any new maps. Need to
%calculate middle of events. Put that in "dates" vector.
m = 1; %index in events vector
iter = 1;   %num. times through outer loop (index of current anomaly map)
while(m < length(event_times_case))% should be length(event_ids)) for a full run
    tic
    n = 1;  %counting vector for num. 3-hrly reports in each event.
    
    %%%%%%%%%%%%%%%%%Take the map closest to the midtime of the event. Averaging is prob a bad idea.
    %Count how many more records are in this event.
    if m < length(event_ids)    %boundary case for end of data
        while(event_ids(m+1) == event_ids(m)) %While the next record is still in this event:
            m = m + 1;
            n = n + 1;
        end
    end
    dates(iter) = event_times_case(ceil(m-(n)/2)); %track times corresponding with maps. This is the median-time of the event.
    
    iter = iter + 1;
    m = m + 1;
    toc
end

%Chop off the first three years of events from the output of
%FZRA_EventTimes:
event_ids = event_ids(event_times >= datetime(1979,1,1) & event_times < datetime(2015,1,1));
event_pct_lightFZRA = event_pct_lightFZRA(event_times >= datetime(1979,1,1) & event_times < datetime(2015,1,1));
event_precip = event_precip(event_times >= datetime(1979,1,1) & event_times < datetime(2015,1,1));
event_precip_std = event_precip_std(event_times >= datetime(1979,1,1) & event_times < datetime(2015,1,1));
event_spd = event_spd(event_times >= datetime(1979,1,1) & event_times < datetime(2015,1,1));
event_spd_std = event_spd_std(event_times >= datetime(1979,1,1) & event_times < datetime(2015,1,1));
event_stationcounts = event_stationcounts(event_times >= datetime(1979,1,1) & event_times < datetime(2015,1,1));

%Create a list of URLs for mean sea level pressure (Pa) and any other variables we'll include:
url_prmsl = [repmat('http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/NARR/monolevel/prmsl.',length(event_ids),1),datestr(event_times_case,'yyyy'),repmat('.nc',length(event_ids),1)];
%Air temp at 2m:
url_airtemp = [repmat('http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/NARR/monolevel/air.2m.',length(event_ids),1),datestr(event_times_case,'yyyy'),repmat('.nc',length(event_ids),1)];
%Geopotential height (m):
url_hgt = [repmat('http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/NARR/pressure/hgt.',length(event_ids),1),datestr(event_times_case,'yyyy'),datestr(event_times_case,'mm'),repmat('.nc',length(event_ids),1)];


%Read in the monolevel stuff:
prmsl_timez = nc_varget(url_prmsl(1,:),'time'); %Figure out new index.
prmsl_timez = datetime(prmsl_timez*3600,'ConvertFrom','epochtime','Epoch','1800-01-01');

%Make anomaly maps for each event:
m = 1; %index in events vector
iter = 1;   %num. times through outer loop (index of current anomaly map)
while(m < length(event_times_case))% should be length(event_ids)) for a full run
    tic
    n = 1;  %counting vector for num. 3-hrly reports in each event.
    
    %%%%%%%%%%%%%%%This was for if we want to make an "event average map" for each event. Bad idea, I think. Just take the map at the midtime of each event.
    %     prmsl_event = [];
    %     %Download a map for each 3-hrs in the event:
    %
    %     %Initial map for first record of event:
    %     indextoplot = find(event_times_case(m) == timez) - 1; %minus one since netCDF uses zero indexing
    %     prmsl_event(:,:,n) = nc_varget(url_prmsl(n,:),'prmsl',[indextoplot yinds xindw],[1 yindn xinde]); %format is (time,y,x)
    %     while(event_ids(m+1) == event_ids(m))
    %         indextoplot = find(event_times_case(m+1) == timez) - 1; %minus one since netCDF uses zero indexing
    %         tic
    %         prmsl_event(:,:,n+1) = nc_varget(url_prmsl(n,:),'prmsl',[0 0 indextoplot],[inf inf 1]); %format is (x,y,time)
    %         toc
    %         m = m + 1;
    %         n = n + 1;
    %     end
    
    %%%%%%%%%%%%%%%%%Take the map closest to the midtime of the event. Averaging is prob a bad idea.
    %Count how many more records are in this event.
    if m < length(event_ids)    %boundary case for end of data
        while(event_ids(m+1) == event_ids(m)) %While the next record is still in this event:
            m = m + 1;
            n = n + 1;
        end
    end
    
    %Now find the index in the prmsl files to download.
    %We'll use the map closest to the median of the event times.
    dates(iter) = event_times_case(ceil(m-(n)/2)); %track times corresponding with maps. This is the median-time of the event.
    
    %Download the times vector for the year in which the midpoint time is
    %and then match it to our timestamp for the index to plot.
    url_prmsl = ['http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/NARR/monolevel/prmsl.',datestr(dates(iter),'yyyy'),'.nc'];
    url_airtemp = ['http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/NARR/monolevel/air.2m.',datestr(dates(iter),'yyyy'),'.nc'];
    url_hgt = ['http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/NARR/pressure/hgt.',datestr(dates(iter),'yyyy'),datestr(dates(iter),'mm'),'.nc'];  
    
    prmsl_timez = nc_varget(url_prmsl,'time'); %Figure out new index.
    prmsl_timez = datetime(prmsl_timez*3600,'ConvertFrom','epochtime','Epoch','1800-01-01');
    hgt_timez = nc_varget(url_hgt,'time'); %Figure out new index for the monthly pressure.
    hgt_timez = datetime(hgt_timez*3600,'ConvertFrom','epochtime','Epoch','1800-01-01');
    
    prmsl_indextoplot = find(dates(iter) == prmsl_timez) - 1; %minus one since netCDF uses zero indexing
    hgt_indextoplot = find(dates(iter) == hgt_timez) - 1; %minus one since netCDF uses zero indexing
    
    %Download these maps to be representative of this event.
    prmsl_anom(:,:,iter) = nc_varget(url_prmsl,'prmsl',[prmsl_indextoplot yinds xindw],[1 yindn xinde]); %format is (time,y,x)
    airtemp(:,:,iter) = nc_varget(url_airtemp,'air',[prmsl_indextoplot yinds xindw],[1 yindn xinde]); %format is (time,y,x)
    hgt850_anom(:,:,iter) = nc_varget(url_hgt850,'hgt',[hgt_indextoplot 6 yinds xindw],[1 1 yindn xinde]); %format is (time,level,y,x)
    hgt1000_anom(:,:,iter) = nc_varget(url_hgt,'hgt',[hgt_indextoplot 0 yinds xindw],[1 1 yindn xinde]); %format is (time,level,y,x)
    hgt500_anom(:,:,iter) = nc_varget(url_hgt,'hgt',[hgt_indextoplot 16 yinds xindw],[1 1 yindn xinde]); %format is (time,level,y,x)
    
    %Subtract out the mean for this month:
    prmsl_anom(:,:,iter) = prmsl_anom(:,:,iter) - prmsl_monthly_avg(:,:,str2double(datestr(dates(iter),'mm')));
    hgt850_anom(:,:,iter) = hgt850_anom(:,:,iter) - hgt850_monthly_avg(:,:,str2double(datestr(dates(iter),'mm')));
    hgt1000_anom(:,:,iter) = hgt1000_anom(:,:,iter) - hgt1000_monthly_avg(:,:,str2double(datestr(dates(iter),'mm')));
    hgt500_anom(:,:,iter) = hgt500_anom(:,:,iter) - hgt500_monthly_avg(:,:,str2double(datestr(dates(iter),'mm')));
    
    
    iter = iter + 1;
    m = m + 1;
    toc
end

%So now we've got our anomaly maps!
prmsl_anom_mb = prmsl_anom*0.01;

%Subtract the heights to get a 1000 - 500 hPa height:
hgt1000500_anom = hgt500_anom - hgt1000_anom;

%load prmsl_anom_19791996
%load maps_prmsl_temp_hgt850_19791996
%load maps_prmsl_temp_hgt1000500_19972014
%load maps_prmsl_temp_hgt1000500_19972014

% %Optional: combine the two periods to make the whole study period for a
% %combined PCA and compare the results with the two halves.
% %Start by saving the current period as a new variable for each:
% hgt1000500_anom1 = hgt1000500_anom;
% prmsl_anom_mb1 = prmsl_anom_mb;
% dates1 = dates;
% airtemp1 = airtemp;
% %Now load the other half of the study period and combine the two:
% hgt1000500_anom = cat(3, hgt1000500_anom1,  hgt1000500_anom);
% prmsl_anom_mb = cat(3, prmsl_anom_mb1, prmsl_anom_mb);
% dates = cat(2, dates1, dates);
% airtemp = cat(3, airtemp1, airtemp);
% %Great! Now save them for later.

load maps_prmsl_temp_hgt1000500_all


%Plot an example map, just to be sure we're cool.
figure(100)
contourf(lon,lat,prmsl_anom_mb(:,:,100),100,'LineColor','none')
colormap(parula)
grid on
title('MSLP during FZRA ')
ylabel('Latitude (deg)')
c = colorbar;
c.Label.String = 'MSL Pressure (mb)';


%% Now PCA THOSE PUPPERS. Reshape the two maps into one long vector.
%So X is a matrix with rows as diff obs and each row being a vector
%representing first the MSL pressure and then the 850mb geopotential
%height.
%X = reshape(prmsl_anom_mb,iter-1,105*99);
Xprmsl = reshape(prmsl_anom_mb,size(prmsl_anom_mb,1)*size(prmsl_anom_mb,2),size(prmsl_anom_mb,3));
%Xhgt850 = reshape(hgt850_anom,size(hgt850_anom,1)*size(hgt850_anom,2),size(hgt850_anom,3));
Xhgt1000500 = reshape(hgt1000500_anom,size(hgt1000500_anom,1)*size(hgt1000500_anom,2),size(hgt1000500_anom,3));
%X = [Xprmsl;Xhgt850];
X = [Xprmsl;Xhgt1000500];


% %Make it into a matrix!
% for z = 1:length(prmsl_anom_mb(1,1,:))
%     row = 1;
%     mat = prmsl_anom_mb(row,:,z);
%     for row = 2:length(prmsl_anom_mb(:,1,1))
%         mat = [mat prmsl_anom_mb(row,:,z)];
%     end
%     X(z,:) = mat;
% end
% 
% X2 = X;
    
numPCs = 10;

%[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(X,'NumComponents',numPCs);
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(X'); %rows are obs, cols are variables

%X_transformed = X*COEFF;
X_transformed = [SCORE*COEFF']';

figure(1)
plot(1:length(EXPLAINED),cumsum(EXPLAINED),'-o')
xlabel('PC')
ylabel('Cumulative Percentage of Variance Explained')
grid on

%Plot the first several PCs.
%Initialize size of PCmat:
PCsprmsl = zeros(105,99,numPCs);
%PCshgt850 = zeros(105,99,numPCs);
PCshgt1000500 = zeros(105,99,numPCs);

for n = 1:numPCs
    %PCs(:,:,n) = reshape(COEFF(:,n),105,99);
    PCsprmsl(:,:,n) = reshape(X_transformed(1:10395,n),size(prmsl_anom_mb,1),size(prmsl_anom_mb,2));
    %PCshgt850(:,:,n) = reshape(X_transformed(10396:end,n),size(hgt850_anom,1),size(hgt850_anom,2));
    PCshgt1000500(:,:,n) = reshape(X_transformed(10396:end,n),size(hgt1000500_anom,1),size(hgt1000500_anom,2));
end

figure(2)
subplot(2,1,1)
contourf(lon,lat,PCsprmsl(:,:,1),100,'LineColor','none')
%colormap(jet)
grid on
title('PC1 of MSL during FZRA ')
ylabel('Latitude (deg)')
c = colorbar;
c.Label.String = 'MSL Pressure (mb)';

subplot(2,1,2)
contourf(lon,lat,PCsprmsl(:,:,2),100,'LineColor','none')
grid on
title('PC2 of MSL during FZRA ')
ylabel('Latitude (deg)')
c = colorbar;
c.Label.String = 'MSL Pressure (mb)';
xlabel('Longitude (deg E of Prime Meridian)')


%SCORES are the expansion coeffs.
figure(3)
plot(dates,SCORE(:,1))
hold on
plot(dates,SCORE(:,2))
xlabel('Time')
ylabel('Expansion Coefficient')
legend('PC1','PC2')

figure(4)
plot(SCORE(:,1),SCORE(:,2))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%Make good PC plots.
%Map boundaries. Make 'em wide
lat_lim = [30 55]; %deg N
lon_lim = [-105 -65]; %deg W

%Import the shapefile:
states = shaperead('usastatehi','UseGeoCoords',true);
provinces = shaperead('province','UseGeoCoords',true);
mexstates = shaperead('shapefiles/mexstates','UseGeoCoords',true);

figure(5)
%Loop through to do all four plots:
for z = 1:4
    subplot(2,2,z)
    worldmap(lat_lim,lon_lim)
    cptcmap('SVS_tempanomaly', 'mapping', 'scaled','flip',true);
    geoshow(states,'FaceColor',[1 1 1]);
    geoshow(provinces,'FaceColor',[1 1 1]);
    geoshow(mexstates,'FaceColor',[1 1 1]);
    %pcolorm(lat,lon,prmsl,'FaceAlpha',0.8)
    
    %Plot MSL pressure:
    [Ca ha] = contourm(lat,lon,PCsprmsl(:,:,z),'LineWidth',1,'LineColor','k');
    pcolorm(lat,lon,PCsprmsl(:,:,z),'FaceAlpha',0.8)
    
    %clabelm(Ca);
    
    %Add an "L" over the low pressure system:
    prmsl_mb = PCsprmsl(:,:,z);
    prmsl_window = prmsl_mb(lat > lat_lim(1) & lat < lat_lim(2) ...
        & lon > lon_lim(1) & lon < lon_lim(2)); %selects only inside our region
    windowmin = min(prmsl_window);
    lattext = lat(prmsl_mb == windowmin);
    lontext = lon(prmsl_mb == windowmin);
    textindex = find(lattext > lat_lim(1) & lattext < lat_lim(2) ...
        & lontext > lon_lim(1) & lontext < lon_lim(2));
    textm(lattext(textindex),lontext(textindex), 'L', 'FontWeight','bold','FontSize',20,'Color','r')
    
    %Add an "H" over any high pressure system (cutoff at 1020 mb)
    windowmax = max(prmsl_window);
    if windowmax > 20 %mb
        lattext = lat(prmsl_mb == max(prmsl_window));
        lontext = lon(prmsl_mb == max(prmsl_window));
        textindex = find(lattext > lat_lim(1) & lattext < lat_lim(2) ...
            & lontext > lon_lim(1) & lontext < lon_lim(2));
        textm(lattext(textindex),lontext(textindex), 'H', 'FontWeight','bold','FontSize',20,'Color','b')
    end
    
    %caxis([windowmin windowmax]/10) %kPa. Sets color ramp to the range of our region.
    caxis([-20 20]) %kPa. Sets color ramp constant for all four plots.
    
    %Add the geopotential heights as dotted lines:
    %[Ca ha] = contourm(lat,lon,PCshgt850(:,:,z),'LineWidth',2,'LineColor','k','LineStyle',':');
    [Ca ha] = contourm(lat,lon,PCshgt1000500(:,:,z),'LineWidth',2,'LineColor','k','LineStyle',':');
    clabelm(Ca);
    
    pause(0.1)
    framem; gridm; tightmap;
    pause(0.1)
    framem; gridm; tightmap;
    pause(0.1)
end

%Add one big colorbar to the whole figure. Just add a small one and adjust
%it.
hp4 = get(subplot(2,2,4),'Position')
c = colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.1  hp4(2)+hp4(3)*2.1])
c.FontSize = 18
c.FontName = 'Gill Sans MT'
xlabel(c,'MSL Pressure Anomaly (mb)')






%% GET THAT K-MEANS:
current_event_def = '43';

load(strcat('maps_',current_event_def,'.mat'))


numclusters = 3;
Xprmsl = reshape(prmsl_anom_mb,size(prmsl_anom_mb,1)*size(prmsl_anom_mb,2),size(prmsl_anom_mb,3));

Xhgt850 = reshape(hgt850_anom,size(hgt850_anom,1)*size(hgt850_anom,2),size(hgt850_anom,3));
X = [Xprmsl;Xhgt850];
% Xhgt1000500 = reshape(hgt1000500_anom,size(hgt1000500_anom,1)*size(hgt1000500_anom,2),size(hgt1000500_anom,3));
% X = [Xprmsl;Xhgt1000500];
[IDX centroids] = kmeans(X', numclusters,'Replicates',50);
%Try it with k-medoids! (Prob will be too slow):
% tic
% [IDX centroids] = kmedoids(X', numclusters);
% toc

% %%%%%%%OR create a self-organizing map with nctool and then do a k-means on
% %%%%%%%THAT. See if we get better performance. If we do, we can redo the
% %%%%%%%station clustering too and bootstrap the kmeans for robustness.
% %nctool
% load test_net_SOM
% X = cell2mat(net.IW)';
% [IDX centroids] = kmeans(X',numclusters);
% %%%%%End of that little bit. Proceed as you were:

%Reshape the vectors into maps:
for m = 1:numclusters
    clustermaps_prmsl(:,:,m) = reshape(centroids(m,1:10395),size(prmsl_anom_mb,1),size(prmsl_anom_mb,2));
    %clustermaps_hgt(:,:,m) = reshape(centroids(m,10396:end),size(hgt850_anom,1),size(hgt850_anom,2));
    clustermaps_hgt(:,:,m) = reshape(centroids(m,10396:end),size(hgt1000500_anom,1),size(hgt1000500_anom,2));
    m
end

% %Plot a silhouette plot to inspect whether or not 3 clusters was best:
% figure(52)
% silhouette(X',IDX)


%Create centroids for the surface air temperature so that we can plot the 0
%deg isotherm on the centroid maps!
for centroidnum = 1:64
    clustermaps_airtempcentroids(:,:,centroidnum) = mean(airtemp(:,:,IDX == centroidnum),3);
end

%%%%%%%%%%%%%%%%%%Good plot:
%Import the shapefile:
states = shaperead('usastatehi','UseGeoCoords',true);
provinces = shaperead('shapefiles/province','UseGeoCoords',true);
mexstates = shaperead('shapefiles/mexstates','UseGeoCoords',true);
lat_lim = [30 55]; %deg N
lon_lim = [-105 -65]; %deg W

figure(6)
x0=2;
y0=2;
width=24;
height=8;
set(gcf,'units','inches','position',[x0,y0,width,height])
dim = [.1 .8 .1 .1];
%annotation('textbox',dim,'String',current_event_def,'FitBoxToText','on','FontSize',48);
%Loop through to do all three cluster plots:
for z = 1:3
    
    subplot(1,3,z)
    worldmap(lat_lim,lon_lim)
    setm(gca,'mapprojection','eqaconic')
    cptcmap('SVS_tempanomaly', 'mapping', 'scaled','flip',true);

    geoshow(states,'FaceColor',[1 1 1]);
    geoshow(provinces,'FaceColor',[1 1 1]);
    geoshow(mexstates,'FaceColor',[1 1 1]);
    %pcolorm(lat,lon,prmsl,'FaceAlpha',0.8)
    
    %Plot MSL pressure:
    [Caz haz] = contourm(lat,lon,clustermaps_prmsl(:,:,z),'LineWidth',1,'LineColor','k');
    pcolorm(lat,lon,clustermaps_prmsl(:,:,z),'FaceAlpha',0.8)
    
    %clabelm(Caz);
    
    %Add an "L" over the low pressure system:
    prmsl_mb = clustermaps_prmsl(:,:,z);
    prmsl_window = prmsl_mb(lat > lat_lim(1) & lat < lat_lim(2) ...
        & lon > lon_lim(1) & lon < lon_lim(2)); %selects only inside our region
    windowmin = min(prmsl_window);
    lattext = lat(prmsl_mb == windowmin);
    lontext = lon(prmsl_mb == windowmin);
    textindex = find(lattext > lat_lim(1) & lattext < lat_lim(2) ...
        & lontext > lon_lim(1) & lontext < lon_lim(2));
    textm(lattext(textindex),lontext(textindex), 'L', 'FontWeight','bold','FontSize',20,'Color','r')
    
    %Add an "H" over any high pressure system (cutoff at 1020 mb)
    windowmax = max(prmsl_window);
    if windowmax > 10 %mb
        lattext = lat(prmsl_mb == max(prmsl_window));
        lontext = lon(prmsl_mb == max(prmsl_window));
        textindex = find(lattext > lat_lim(1) & lattext < lat_lim(2) ...
            & lontext > lon_lim(1) & lontext < lon_lim(2));
        textm(lattext(textindex),lontext(textindex), 'H', 'FontWeight','bold','FontSize',20,'Color','b')
    end
    
    caxis([-20 20]) %kPa. Sets color ramp to the range of our region.
    

    pause(1)
    framem; gridm; tightmap;
    pause(1)
    pause(1)
    framem; gridm; tightmap;
    pause(1)
    
    %Add the upper air geopotential heights as dotted lines:
    
    [Ca ha] = contourm(lat,lon,clustermaps_hgt(:,:,z),'LineWidth',2,'LineColor','k','LineStyle',':');
    clabelm(Ca);
    
    %Lastly, add a 0 deg isotherm in a nice bold black:
    [Cb hb] = contourm(lat,lon,clustermaps_airtempcentroids(:,:,z),[273.15 273.15],'LineWidth',2,'LineColor','k');
    
end

%Add one big colorbar to the whole figure. Just add a small one and adjust
%it.
hp4 = get(subplot(1,3,3),'Position')
c = colorbar('Position', [hp4(1)+hp4(3)+0.012  hp4(2)*1.5  0.025  hp4(2)+hp4(3)*2.5])
c.FontSize = 18
c.FontName = 'Helvetica'
xlabel(c,'MSL Pressure Anomaly (mb)')



% %% Add new clusters by computing distance to the centroids. Then analyze the change in prevalence over time. 
% % load prmsl_anom_19791996
% % load IDX_and_X
% %load data19791996
% %load prmsl_anom_19982014 %has _recent appended to all variable names
% %load data19972014  %has _recent appended to all variable names
% 
% %Loop through each new event, and assign it to the smallest distance
% %centroid:
% for m = 1:size(X_recent,2)
%     for centroidnum = 1:3
%         %distance(centroidnum) = sum(sqrt(X_recent(:,m).^2 + centroids(centroidnum,:)'.^2));
%         distance(centroidnum) = norm(X_recent(:,m) - centroids(centroidnum,:)');
%     end
%     [Y clusterpick] = min(distance);
%     IDX_recent_orig(m) = clusterpick;   %vector of classifications into the original cluster from the 1979-1996 period
% end
% 
% %Plot this stuff
% figure(1000)
% scatter(dates,IDX)
% hold on
% grid on
% scatter(dates_recent,IDX_recent_orig)
% xlabel('Time (years)')
% ylabel('Cluster')
% 
% figure(1001)
% % that = [sum(IDX == 1)/586*100  sum(IDX_recent_orig == 1)/625*100; ...
% %     sum(IDX == 2)/586*100  sum(IDX_recent_orig == 2)/625*100; ...
% %     sum(IDX == 3)/586*100  sum(IDX_recent_orig == 3)/625*100]
% that = [sum(IDX == 1)  sum(IDX_recent_orig == 1); ...
%     sum(IDX == 2)  sum(IDX_recent_orig == 2); ...
%     sum(IDX == 3)  sum(IDX_recent_orig == 3)]
% 
% bar(that)
% ylabel('Number of Events')
% xlabel('Cluster')
% legend('1979-1996','1997-2014')
% grid on
% 
% 
% 
% 
% 








%% MAP PLOTTING OF CLUSTER STUFF
load a_all.mat
load b.mat %contains b, the full dataset of FZRA events for each station. *fixed b file. coords are fixed too. can prob delete b_fixedcoords
stationnames = fieldnames(b);
% %Specify a map area of interest (the bounds of our map):
% lat_lim = [38 50]; %deg N
% lon_lim = [-95 -71]; %deg W



% %% Plot pie charts at each location to show the relative number of storms in each of the three bins at each location.
% %First make clusterpct, a 97x3 matrix that shows the percentage of synoptic
% %forcing categories that makes up the storms at each station:
% %Loop through all events, matching their dates up with output from
% %FZRA_EventTimes and then giving a count to each station that participated
% %in the event.
% clusterpct = zeros(length(stationnames),numclusters);
% 
% %This only assumes that dates are a subset of event_times
% for m = 1:length(dates)
%     dateindex = event_ids(dates(m) == event_times);   %finds the matching event ID
%     %Now add a count representing one event participated in for each
%     %station that participated in this event, under the right category:
%     clusterpct(~~event_stationcounts(dateindex,:),IDX(m)) = clusterpct(~~event_stationcounts(dateindex,:),IDX(m)) + 1;
% end
% 
% %load clusters_n_k5stationclusters
% %load clusters_n_k5stationclusters_BETTER
% load clusters_n_k5stationclusters_BETTER_plus_NYC_6thcluster
% 
% %We could normalize them all by percentage, but we could just leave them
% %and use their sum as a scaling for the size of the pie charts if we wanna
% %be obnoxious.
% 
% %Plot the points in 3d so we can see if the clusters make sense:
% % figure(20)
% % plot3(clusterpct(:,1),clusterpct(:,2),clusterpct(:,3),'o','LineStyle','none')
% % xlabel('Arctic High & Cold Air Damming')
% % ylabel('Cyclone/Anticyclone')
% % zlabel('Occluded Front & Cold Air Trapping')
% % axis vis3d
% % view(90,0)
% % grid on
% 
% figure(100)
% worldmap(lat_lim,lon_lim)
% colormap(parula)
% %caxis([-40, 40])
% geoshow(states,'FaceColor',[1 1 1])
% hold on
% geoshow(provinces,'FaceColor',[1 1 1])
% framem('ffacecolor',[.5 .7 .9]); %shows water as blue
% %plotm(a.StationLocations(:,1),a.StationLocations(:,2),'ko') %plots measurement stations
% 
% %Plot coloring stations by how many storms they've participated in:
% %scatterm(a.StationLocations(:,1),a.StationLocations(:,2),250,sum(event_stationcounts),'filled')
% 
% %colorbar
% hold on
% title('Category of Synoptic Storm Type by Location')
% 
% %Loop through each station, creating and placing a pie chart:
% for station = 1:length(stationnames)
%     p = pie(clusterpct(station,:));     %where clusterpct is a 97 x 3 matrix
%     %p(1).Vertices
%     sc = 0.45; % <scaling factor (different for each city)
%     % Loop through each slice of the pie:
%     for k = 1:2:length(p)
%         % x,y coordinates of a slice:
%         tmp = p(k).Vertices;
%         % Scale the size of the slice:
%         tmp = tmp*sc;
%         % Place the center of the pie on its station:
%         tmp(:,1) = tmp(:,1)+b.(stationnames{station}).coords(1);
%         tmp(:,2) = tmp(:,2)+b.(stationnames{station}).coords(2);
%         
%         switch k
%             case 1 
%                 patchm(tmp(:,1),tmp(:,2),[1 0.1 0.1]);
%             case 3
%                 patchm(tmp(:,1),tmp(:,2),[0.75 0.75 0.1]);
%             case 5
%                 patchm(tmp(:,1),tmp(:,2),[0.1 0.1 0.8]);
%         end
%     end
%     
% end



%% Now, cluster the stations based off their breakdowns. Will they give regional identities?
%Prob not very well. But ya gotta try.
numclusters = 5;
tic
[IDX_stations centroids_stations] = kmeans(clusterpct, numclusters,'Replicates',1000);
toc

% %%%%%%%%%%%%%OR we can do it with a time bias where we let time be included
% %%%%%%%%%%%%%as a factor in the clustering. So not only will the stations
% %%%%%%%%%%%%%be clustered by the percentage of participation overall in
% %%%%%%%%%%%%%the three storm types, but they'll be selected to be near
% %%%%%%%%%%%%%stations that they experienced each storm with. Wait, this is
% %%%%%%%%%%%%%a whole different way to classify them. It's just based off of
% %%%%%%%%%%%%%which stations had FZRA at the same time. Whatever, let's try
% %%%%%%%%%%%%%it:
% numclusters = 6;
% [IDX_stations centroids_stations] = kmeans(event_stationcounts', numclusters);
% clusterpct = event_stationcounts';

% figure(49)
% silhouette(event_stationcounts',IDX_stations)

%eva = evalclusters(event_stationcounts,'kmeans','CalinskiHarabasz','KList',[1:6])

%Plot a silhouette plot to inspect whether or not 5 clusters was best:
figure(51)
silhouette(clusterpct,IDX_stations)

% %Now we're gonna create our own SIXTH CATEGORY for the NYC-Long Island
% %area:
% cat6indices = [24 26 27 97]; %these are KLGA, KHPN, KISP, and KJFK (in that order)
% IDX_stations(cat6indices) = 6;
% %Add a nice new color to colorz for plotting.
% colorz(6,:) = [255,127,80]/255;

%Plot the points in 3d so we can see if the clusters make sense:
load colorz
figure(23)
for k = 1:numclusters + 1  %plus one for our additional NYC Category
    %colorz(k,:) = rand(1,3); %Use the loaded colors if you've loaded results
    %plot3(clusterpct(IDX_stations==k,1),clusterpct(IDX_stations==k,2),clusterpct(IDX_stations==k,3),'o','color',colorz(k,:),'LineStyle','none','LineWidth',3)
    %plot(clusterpct(IDX_stations==k,3),clusterpct(IDX_stations==k,1),'o','color',colorz(k,:),'LineStyle','none','LineWidth',3)
    subplot(1,3,1)
    plot(clusterpct(IDX_stations==k,1),clusterpct(IDX_stations==k,3),'o','color',colorz(k,:),'LineStyle','none','LineWidth',3)
    hold on
    xlabel('1')
    ylabel('3')
%     xlabel('Arctic High & Cold Air Damming')
%     ylabel('Cyclone/Anticyclone')
%     zlabel('Occluded Front & Cold Air Trapping')
    xlim([0 90])
    %axis vis3d
    grid on
    
    subplot(1,3,2)
    plot(clusterpct(IDX_stations==k,2),clusterpct(IDX_stations==k,1),'o','color',colorz(k,:),'LineStyle','none','LineWidth',3)
    hold on
    xlabel('2')
    ylabel('1')
%     xlabel('Arctic High & Cold Air Damming')
%     ylabel('Cyclone/Anticyclone')
%     zlabel('Occluded Front & Cold Air Trapping')
    xlim([0 90])
    grid on
    
    subplot(1,3,3)
    plot(clusterpct(IDX_stations==k,3),clusterpct(IDX_stations==k,2),'o','color',colorz(k,:),'LineStyle','none','LineWidth',3)
    hold on
    xlabel('3')
    ylabel('2')
%     xlabel('Arctic High & Cold Air Damming')
%     ylabel('Cyclone/Anticyclone')
%     zlabel('Occluded Front & Cold Air Trapping')
    xlim([0 90])
    grid on
    
    
end

%plot the stations by their cluster identities:
lat_lim = [37.5 50]; %deg N
lon_lim = [-95 -71]; %deg W
figure(2)
worldmap(lat_lim,lon_lim) %plots empty axes
%title('Regional Classification from Synoptic Clustering');
geoshow(states,'FaceColor',[1 1 1])
geoshow(provinces,'FaceColor',[1 1 1])
framem('ffacecolor',[.5 .7 .9]); %shows water as blue
for k = 1:numclusters + 1  %plus one for our additional NYC Category
    plotm(a.StationLocations(IDX_stations==k,1),a.StationLocations(IDX_stations==k,2),'+','color',colorz(k,:),'MarkerSize',10,'LineWidth',4)
    hold on
end


