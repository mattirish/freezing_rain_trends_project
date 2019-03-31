%% FZRADurationAnalysis.m
% Analyzes spatial and temporal trends in the duration of FZRA events
% throughout the Greate Lakes
% Matt Irish
% March 2017

% Uses the 'b' struct.
clear

load a_all.mat
load b.mat

%% Loop through all stations, years, and months mayb, adding up duration counts.

stationnames = fieldnames(b);

% Initialize EZ access 97*39 mean and median event duration matrices:
meandurations = zeros(97,39);
mediandurations = meandurations;

%Outer loop: loop through stations.
for m = 1:length(stationnames)
    stn = b.(stationnames{m})       %sets current station
    
    %Inner loop: Go through each event.
    eventdur = 1;       %initialize event duration counter
    durationrecords = [];
    durationyr = [];
    for n = 1:length(stn.YR)-1
        if stn.HR(n+1) == stn.HR(n) + 1 ...         %next hour is the next record
                & stn.DA(n+1) == stn.DA(n)         %next record is the same day
            eventdur = eventdur + 1;                %add an hour to the current event
        elseif stn.HR == 23 ...                     %hour is 23:00 (special case)
                & stn.HR(n+1) == 0 ...             %midnight is next record
                & stn.DA(n+1) == stn.DA(n) + 1     %next record is next day
            eventdur = eventdur + 1;                %add an hour to the current event
        else
            durationrecords(end+1) =  eventdur;      %event is over. reset counter and add
            durationyr(end+1) = stn.YR(n);          %record year this event took place in
            eventdur = 1;                           %reset event duration counter
        end
    end
    
    %Compute the mean and median duration for this station, and make a
    %matrix of event durations for each year:
    for yr = 1976:2014
        meandurations(m,yr-1975) = mean(durationrecords(durationyr == yr));
        mediandurations(m,yr-1975) = median(durationrecords(durationyr == yr));
        
        %Annual event duration distribution. 1 to 9 hours and then 10th
        %element is events lasting more than 10 hours:
        for hr = 1:9
            b.(stationnames{m}).durationdistrib(yr-1975,hr) = sum(durationyr == yr & durationrecords == hr);
        end
        b.(stationnames{m}).durationdistrib(yr-1975,10) = sum(durationyr == yr & durationrecords >= 9);  %>9 hours
        
        %Make each year's distribution into a percent of total records for that year:
        b.(stationnames{m}).durationdistrib(yr-1975,:) = b.(stationnames{m}).durationdistrib(yr-1975,:)./sum(durationyr == yr)*100; 
    end
    
    %Replicate the Cortinas 2000 bar plot and then do an updated one.
    %To do this, we need station average distribs for each period:
    %nanszeroed = b.(stationnames{m}).durationdistrib;
    %nanszeroed(isnannanszeroed) = 0;
    durationdistrib76to90(m,:) = sum(b.(stationnames{m}).durationdistrib(1:15,:),'omitnan')/15;
    durationdistrib00to14(m,:) = sum(b.(stationnames{m}).durationdistrib(25:39,:),'omitnan')/15;
end

%Average all the stations to replicate the Cortinas 2000 bar plot and the
%last 15 years updated plot:
durationdistrib76to90 = mean(durationdistrib76to90);
durationdistrib00to14 = mean(durationdistrib00to14);

%Save duration stuff:
%save('durationanalysis.mat','meandurations','mediandurations','durationdistrib76to90','durationdistrib00to14')


%% Good job Matt! Now let's actually make these plots. 
%% Bar plots for Cortinas replicate and update:
durations = [durationdistrib76to90;durationdistrib00to14]';
figure(1)
bar(durations)
Labels = {'1','2','3','4','5','6','7','8','9','>9'};
set(gca, 'XTick', 1:10, 'XTickLabel', Labels);
legend('1976-1990','2000-2014')
grid on
title('Average Duration of FZRA Events throughout Region')
xlabel('Number of Consecutive Hourly Reports')
ylabel('Percentage of FZRA Events')

%% Map of trends in average event duration. Will we see a decrease in the South? Increase in the North?
lat_lim = [39.5 50]; %deg N
lon_lim = [-94 -72.5]; %deg W

%Import shapefiles:
states = shaperead('usastatehi','UseGeoCoords',true);
provinces = shaperead('province','UseGeoCoords',true);

t = 1976:2014;
thresh = 0.1;       %set statistical signifance threshold

%Run statistical analysis for each station:
for m = 1:97
    [taub(m) tau(m) h(m) sig(m) Z(m) S(m) sigma(m) sen(m) n senplot cilower(m) ciupper(m)] = ktaub([t;meandurations(m,t-1975)]',thresh,1);
end

%Decadal trend in percent relative to 1976-1985 baseline:
sen_percent = sen'./nanmean(meandurations(:,1:10),2)*100*10;
 
figure(2)
worldmap(lat_lim,lon_lim)

%cptcmap('SVS_tempanomaly', 'mapping', 'scaled');
c = colorbar; 
caxis([-15, 15])

geoshow(states,'FaceColor',[.75 .75 .75])
hold on
geoshow(provinces,'FaceColor',[.75 .75 .75])
framem('ffacecolor',[.5 .7 .9]); %shows water as blue
%plotm(a.StationLocations(:,1),a.StationLocations(:,2),'r+') %plots measurement stations

sizecoeff = 700;
thresh = 1;     %set significance threshold for display
scatterm(a.StationLocations(sig < thresh,1),a.StationLocations(sig < thresh,2),500*(1-sig(sig < thresh)),sen_percent(sig < thresh),'filled')
c = colorbar
c.FontSize = 16
ylabel = get(c,'YTickLabel');
perc = repmat('%',size(ylabel,1),1);
ylabel = strcat(ylabel,cellstr(perc));
set(c,'YTickLabel',ylabel);
xlabel(c,'Change in Mean Event Duration per Decade')
hold on
title('Trends in Mean FZRA Event Duration Rel. to 1976-85 Baseline','FontSize',18)




figure(3)
surf(meandurations)


figure(4)
plot(sig)
hold on
plot(0.05*ones(size(sig)))
sum(sig < 0.05)/numel(sig)