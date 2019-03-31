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
% Matt Irish 2017

clear 

load monthlies_correctedfordomainsize %archive of MSL maps for all 456 months and 12 averages, one for each month!


%% Load the system boundaries.
%Load the lat, lon, x, and y, (along with cfrzr, from 1979), the subregion
%domain I chose. Routine is at the beginning of
%FZRA_SynopticWeatherTyping.m.
load latslons_subregion_for_PCA


%% Reduce the domain of the monthlies down to our chosen region.
% prmsl_monthly_avg = prmsl_monthly_avg(xindw+1:xindw+xinde,yinds+1:yinds+yindn,:);
% %Rearrange so that it goes y,x,time:
% prmsl_monthly_avg = permute(prmsl_monthly_avg,[2,1,3]);
% %Okay, I think we're good to go now. The dim order is y by x by time.

%% Download surface MSL maps for all the events of 2014 and create anomaly 
% maps of each by subtracting the event average from the all-time monthly
% average.

%Call FZRA_EventTimes to give us a 3-hourly NARR-ready log of all events:
timestep = 3;               %rounds to every three hours
min_reports_per_event = 3;  %min no. reports that constitute an event
max_nonevent_hrs = 3;       %allow up to 4 hours between events

%[event_times event_ids] = FZRA_EventTimes(timestep, min_reports_per_event, max_nonevent_hrs);
load allevents_ids_&_times

%Select only dates for the decided time window:
event_ids = event_ids(event_times >= datetime(1997,1,1) & event_times < datetime(2015,1,1));

event_times_case = event_times(event_times >= datetime(1997,1,1) & event_times < datetime(2015,1,1));
event_str = datestr(event_times_case);

%Create a URL for mean sea level pressure (Pa), just to download the time variable:
url_prmsl = ['http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/NARR/monolevel/prmsl.',datestr(event_times_case(1),'yyyy'),'.nc'];

%Read in the monolevel stuff:
timez = nc_varget(url_prmsl,'time'); %Figure out new index.
timez = datetime(timez*3600,'ConvertFrom','epochtime','Epoch','1800-01-01');

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
    prmsl_anom_dates(iter) = event_times_case(ceil(m-(n)/2)); %track times corresponding with maps. This is the median-time of the event.
    
    %Download the times vector for the year in which the midpoint time is
    %and then match it to our timestamp for the index to plot.
    url_prmsl = ['http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/NARR/monolevel/prmsl.',datestr(prmsl_anom_dates(iter),'yyyy'),'.nc'];
    %Air temp at 2m:
    url_airtemp = ['http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/NARR/monolevel/air.2m.',datestr(prmsl_anom_dates(iter),'yyyy'),'.nc'];
    %Geopotential height (m):
    url_hgt = ['http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/NARR/pressure/hgt.',datestr(prmsl_anom_dates(iter),'yyyy'),datestr(prmsl_anom_dates(iter),'mm'),'.nc'];
    
    timez = nc_varget(url_prmsl,'time'); %Figure out new index.
    timez = datetime(timez*3600,'ConvertFrom','epochtime','Epoch','1800-01-01');
    
    indextoplot = find(prmsl_anom_dates(iter) == timez) - 1; %minus one since netCDF uses zero indexing
    
    %Download this map to be representative of this event.
    prmsl_anom(:,:,iter) = nc_varget(url_prmsl,'prmsl',[indextoplot yinds xindw],[1 yindn xinde]); %format is (time,y,x)
    hgt500_anom(:,:,iter) = nc_varget(url_hgt,'hgt',[indextoplot 16 yinds xindw],[1 1 yindn xinde]); %format is (time,level,y,x). 17th level is 500mb.
    airtemp(:,:,iter) = nc_varget(url_airtemp,'air',[indextoplot yinds xindw],[1 yindn xinde]); %format is (time,y,x)
    
    %Subtract out the mean for this month:
    prmsl_anom(:,:,iter) = prmsl_anom(:,:,iter) - prmsl_monthly_avg(:,:,str2double(datestr(prmsl_anom_dates(iter),'mm')));
    hgt500_anom(:,:,iter) = hgt500_anom(:,:,iter) - hgt500_monthly_avg(:,:,str2double(datestr(prmsl_anom_dates(iter),'mm')));
    
    
    iter = iter + 1;
    m = m + 1;
    toc
end

%So now we've got our anomaly maps!
prmsl_anom_mb = prmsl_anom*0.01;

%load prmsl_anom_19791996


%Plot an example map, just to be sure we're cool.
figure(100)
contourf(lon,lat,prmsl_anom_mb(:,:,100),100,'LineColor','none')
colormap(parula)
grid on
title('MSLP during FZRA ')
ylabel('Latitude (deg)')
c = colorbar;
c.Label.String = 'MSL Pressure (mb)';


%% Now PCA THOSE PUPPERS:
%X = reshape(prmsl_anom_mb,iter-1,105*99);
X = reshape(prmsl_anom_mb,size(prmsl_anom_mb,1)*size(prmsl_anom_mb,2),size(prmsl_anom_mb,3));

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
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(X);

prmsl_anom_mb_transformed = X*COEFF;

figure(1)
plot(1:iter-1,cumsum(EXPLAINED))
xlabel('PC')
ylabel('Cumulative Percentage of Variance Explained')
grid on

%Plot the first two PCs.
%Initialize size of PCmat:
PCs = zeros(105,99,numPCs);
for n = 1:numPCs
    %PCs(:,:,n) = reshape(COEFF(:,n),105,99);
    PCs(:,:,n) = reshape(prmsl_anom_mb_transformed(:,n),size(prmsl_anom_mb,1),size(prmsl_anom_mb,2));
end
%Now PCs is a 31x180x10 double. It's all the PCs maps.

figure(2)
subplot(2,1,1)
contourf(lon,lat,PCs(:,:,1),100,'LineColor','none')
%colormap(jet)
grid on
title('PC1 of MSL during FZRA ')
ylabel('Latitude (deg)')
c = colorbar;
c.Label.String = 'MSL Pressure (mb)';

subplot(2,1,2)
contourf(lon,lat,PCs(:,:,2),100,'LineColor','none')
grid on
title('PC2 of MSL during FZRA ')
ylabel('Latitude (deg)')
c = colorbar;
c.Label.String = 'MSL Pressure (mb)';
xlabel('Longitude (deg E of Prime Meridian)')


%SCORES are the expansion coeffs.
figure(3)
plot(prmsl_anom_dates,SCORE(:,1))
hold on
plot(prmsl_anom_dates,SCORE(:,2))
xlabel('Time')
ylabel('Expansion Coefficient')
legend('PC1','PC2')

figure(4)
plot(SCORE(:,1),SCORE(:,2))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%Make good PC plots.
%Map boundaries. Make 'em wide
lat_lim = [25 60]; %deg N
lon_lim = [-105 -60]; %deg W

%Import the shapefile:
states = shaperead('usastatehi','UseGeoCoords',true);
provinces = shaperead('province','UseGeoCoords',true);

figure(5)
%Loop through to do all four plots:
for z = 1:4
    subplot(2,2,z)
    worldmap(lat_lim,lon_lim)
    cptcmap('SVS_tempanomaly', 'mapping', 'scaled','flip',true);
    geoshow(states,'FaceColor',[1 1 1]);
    geoshow(provinces,'FaceColor',[1 1 1]);
    %pcolorm(lat,lon,prmsl,'FaceAlpha',0.8)
    
    %Plot MSL pressure:
    [Ca ha] = contourm(lat,lon,PCs(:,:,z)/10,'LineWidth',1,'LineColor','k');
    pcolorm(lat,lon,PCs(:,:,z)/10,'FaceAlpha',0.8)
    
    clabelm(Ca);
    
    %Add an "L" over the low pressure system:
    prmsl_mb = PCs(:,:,z);
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
xlabel(c,'MSL Pressure (mb)')






%% GET THAT K-MEANS:
numclusters = 3;
[IDX centroids] = kmeans(X', numclusters);

%Reshape the vectors into maps:
for m = 1:numclusters
    clustermaps(:,:,m) = reshape(centroids(m,:),size(prmsl_anom_mb,1),size(prmsl_anom_mb,2));
end

%Plot one of these bad boyz without the mapping toolbox:
figure(6)
contourf(lon,lat,clustermaps(:,:,1),100,'LineColor','none')
grid on
ylabel('Latitude (deg)')
c = colorbar;


%%%%%%%%%%%%%%%%%%Good plot:
figure(6)
%Loop through to do all three cluster plots:
for z = 1:3
    subplot(2,2,z)
    worldmap(lat_lim,lon_lim)
    cptcmap('SVS_tempanomaly', 'mapping', 'scaled');

    geoshow(states,'FaceColor',[1 1 1]);
    geoshow(provinces,'FaceColor',[1 1 1]);
    %pcolorm(lat,lon,prmsl,'FaceAlpha',0.8)
    
    %Plot MSL pressure:
    [Caz haz] = contourm(lat,lon,clustermaps(:,:,z),'LineWidth',1,'LineColor','k');
    pcolorm(lat,lon,clustermaps(:,:,z),'FaceAlpha',0.8)
    
    clabelm(Caz);
    
    %Add an "L" over the low pressure system:
    prmsl_mb = clustermaps(:,:,z);
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
end

%Add one big colorbar to the whole figure. Just add a small one and adjust
%it.
hp4 = get(subplot(2,2,4),'Position')
c = colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.1  hp4(2)+hp4(3)*2.1])
c.FontSize = 18
c.FontName = 'Gill Sans MT'
xlabel(c,'MSL Pressure (mb)')



%% Add new clusters by computing distance to the centroids. Then analyze the change in prevalence over time. 
load prmsl_anom_19791996
load IDX_and_X
load prmsl_anom_19982014 %has _recent appended to all variable names

%Loop through each new event, and assign it to the smallest distance
%centroid:
for m = 1:size(X_recent,2)
    for centroidnum = 1:3
        %distance(centroidnum) = sum(sqrt(X_recent(:,m).^2 + centroids(centroidnum,:)'.^2));
        distance(centroidnum) = norm(X_recent(:,m) - centroids(centroidnum,:)');
    end
    [Y clusterpick] = min(distance);
    IDX_recent_orig(m) = clusterpick;   %vector of classifications into the original cluster from the 1979-1996 period
end

%Plot this stuff
figure(1000)
scatter(prmsl_anom_dates,IDX)
hold on
grid on
scatter(prmsl_anom_dates_recent,IDX_recent_orig)

figure(1001)
that = [sum(IDX == 1)/586*100  sum(IDX_recent_orig == 1)/548*100; ...
    sum(IDX == 2)/586*100  sum(IDX_recent_orig == 2)/548*100; ...
    sum(IDX == 3)/586*100  sum(IDX_recent_orig == 3)/548*100]
bar(that)
ylabel('% of Events')
legend('1979-1997','1998-2014')
grid on




