function [ event_times, event_ids, event_spd, event_spd_std, event_precip, event_precip_std, event_pct_lightFZRA, event_stationcounts, nonevent_times, nonevent_stations ] = FZRA_EventTimes( timestep, min_reports_per_event, max_nonevent_hrs )
%FZRA_EventTimes: Outputs a vector of FZRA event times for a given input date
%range and minimum number of simultaneous reports to define an event. Also
%outputs a struct showing the number of reports at each station. Order of
%stations is the same as that in b.mat.
%   Inputs:
%       -timestep: will return record every timestep hours (1, 3, or 12
%        for actual hourly, NARR, and rawinsonde analysis, respectively)
%       -min_reports_per_event: min threshold of events to consider an "event"
%       -max_nonevent_hrs: max num. hours that can intervene between FZRA
%        reports within one event
%
%   Outputs:
%       -event_times: unique event times in ascending order
%       -event_id: unique event identifier corresponding to each
%        element of event_times
%       -event_spd: mean wind speed during event at stations with reports
%       -event_spd_std: standard dev of mean wind speed during event at stations with reports
%       -event_precip: mean total liquid precip during event at stations
%       with reports (not including any intervening hours)
%       -event_spd_std: standard dev of mean precip during event at stations with reports
%       -event_pct_lightFZRA: percent of reports for storm saying FZRA is light
%       -event_stationcounts: length(end(event_id)) x length(fieldnames(b)) matrix
%       showing the number of reports at each station for each event.
%       (Replaced "nonevents" double with a vector of datetimes (nonevent_times) and a vector of nonevent stations (nonevent_stations) matrix in August 2018:)
%       -nonevent_times: datetimes of events recorded for which an
%       event wasn't triggered (min threshold of reports wasn't met).
%       -nonevent_stations: station numbers corresponding with the nonevent
%       reports in nonevent_times
%
%--------------------------------------------------------------------------
%% Load all times at all stations where FZRA was reported and sort them 
%  into ascending order, keeping duplicates.

%Load historical freezing rain reports
load b.mat %contains b, the full dataset for each station.

%Following part is based on FZRADurationAnalysis.m:
%Load station names:
stationnames = fieldnames(b);

%Initialize total records datetime vector with first station:
clear total_records             %this will be a vector of all reports of FZRA as time stamps
clear total_stationnums_of_records    %for each record, keep the station number (as an index of stationnames, not the actual number)
clear total_spd                %for each record, keep the wind and precip
clear total_precip 
clear total_islightFZRA        %for each record, keep the intensity from both the manual and auto obs (if any). 1 = light, 0 = mod or heavy

%Initialize to create the variables:
m = 1;
total_records = datetime(b.(stationnames{1}).YR,b.(stationnames{1}).MO,...
        b.(stationnames{1}).DA,b.(stationnames{1}).HR,...
        zeros(size(b.(stationnames{1}).HR)),zeros(size(b.(stationnames{1}).HR)));
total_stationnums_of_records = m*ones(length(b.(stationnames{m}).YR),1);
total_spd = b.(stationnames{m}).SPD;
total_precip = b.(stationnames{m}).PCP01;

%For all reports that are deemed light FZRA, give them a value of 1.
%Otherwise, zero.
total_islightFZRA = -999*ones(length(b.(stationnames{m}).MW1),1); %initialize records to "other intensity category" = -999. Hopefully there won't be any of these.
total_islightFZRA(b.(stationnames{m}).MW1 == 66 | b.(stationnames{m}).MW2 == 66 |  ...  %check if any of them are 66, meaning "Rain, freezing, slight" as opposed to 67, "Rain, freezing, moderate or heavy"
   b.(stationnames{m}).MW3 == 66 | b.(stationnames{m}).MW4 == 66 |  ...
   b.(stationnames{m}).AW1 == 64 | b.(stationnames{m}).AW2 == 64 |  ...  %same as before but light is 64 for automated obs.
   b.(stationnames{m}).AW3 == 64 | b.(stationnames{m}).AW4 == 64) = 1;
%Set everything that's moderate or heavy to zero:
total_islightFZRA(b.(stationnames{m}).MW1 == 67 | b.(stationnames{m}).MW2 == 67 |  ...  %67, "Rain, freezing, moderate or heavy"
   b.(stationnames{m}).MW3 == 67 | b.(stationnames{m}).MW4 == 67 |  ...
   b.(stationnames{m}).AW1 == 65 | b.(stationnames{m}).AW2 == 65 |  ...  %65, "Rain, freezing, moderate"
   b.(stationnames{m}).AW3 == 65 | b.(stationnames{m}).AW4 == 65 |  ...
   b.(stationnames{m}).AW1 == 66 | b.(stationnames{m}).AW2 == 66 |  ...  %66, "Rain, freezing, heavy"
   b.(stationnames{m}).AW3 == 66 | b.(stationnames{m}).AW4 == 66) = 0;

    
%Loop through stations to create one big chronological list of
%reports and their parameters:
for m = 2:length(stationnames)
    records = datetime(b.(stationnames{m}).YR,b.(stationnames{m}).MO,...
        b.(stationnames{m}).DA,b.(stationnames{m}).HR,...
        zeros(size(b.(stationnames{m}).HR)),zeros(size(b.(stationnames{m}).HR)));
    stationnums_of_records = m*ones(length(b.(stationnames{m}).YR),1); %Make another vector of the station numbers to say which reports were
    %from which stations. We'll make this into a matrix so we can ascribe
    %reports to individual events and non-events.
    wind = b.(stationnames{m}).SPD;
    precip = b.(stationnames{m}).PCP01;
    
    
    %For all reports that are deemed light FZRA, give them a value of 1.
    %Otherwise, zero.
    lightfzra = -999*ones(length(b.(stationnames{m}).MW1),1); %initialize records to "other intensity category" = -999. Hopefully there won't be any of these.
    lightfzra(b.(stationnames{m}).MW1 == 66 | b.(stationnames{m}).MW2 == 66 |  ...  %check if any of them are 66, meaning "Rain, freezing, slight" as opposed to 67, "Rain, freezing, moderate or heavy"
        b.(stationnames{m}).MW3 == 66 | b.(stationnames{m}).MW4 == 66 |  ...
        b.(stationnames{m}).AW1 == 64 | b.(stationnames{m}).AW2 == 64 |  ...  %same as before but light is 64 for automated obs.
        b.(stationnames{m}).AW3 == 64 | b.(stationnames{m}).AW4 == 64) = 1;
    %Set everything that's moderate or heavy to zero:
    lightfzra(b.(stationnames{m}).MW1 == 67 | b.(stationnames{m}).MW2 == 67 |  ...  %67, "Rain, freezing, moderate or heavy"
        b.(stationnames{m}).MW3 == 67 | b.(stationnames{m}).MW4 == 67 |  ...
        b.(stationnames{m}).AW1 == 65 | b.(stationnames{m}).AW2 == 65 |  ...  %65, "Rain, freezing, moderate"
        b.(stationnames{m}).AW3 == 65 | b.(stationnames{m}).AW4 == 65 |  ...
        b.(stationnames{m}).AW1 == 66 | b.(stationnames{m}).AW2 == 66 |  ...  %66, "Rain, freezing, heavy"
        b.(stationnames{m}).AW3 == 66 | b.(stationnames{m}).AW4 == 66) = 0;
    
    %Concatenate these records unto the region-wide vector:
    total_records = [total_records; records];
    total_stationnums_of_records = [total_stationnums_of_records; stationnums_of_records];
    total_spd = [total_spd; wind];
    total_precip = [total_precip; precip];
    total_islightFZRA = [total_islightFZRA; lightfzra];

    
end

%Now sort 'em all so they're chronological!
[total_records, order] = sort(total_records); %ascending
total_stationnums_of_records = total_stationnums_of_records(order);
total_spd = total_spd(order);
total_precip = total_precip(order);
total_islightFZRA = total_islightFZRA(order);




%--------------------------------------------------------------------------
%% Now apply the minimum simultaneous reports threshold to define the number
%of simultaneous reports needed to classify it as a FZRA "event" and then
%create a companion vector that classifies each event in increasing order.

%Work through each event.
%initialize.
event_times = [];
event_ids = [];
event_spd = [];
event_spd_std = [];
event_precip = [];
event_precip_std = [];
event_pct_lightFZRA = [];
event_stationcounts = [];
nonevent_times = [];
nonevent_stations = [];
nonevent_stationlist = [];

index = 1; %position in total_records
id = 1;     %initialize event ID

%This will loop once per event (or nonevent):
while(index <= length(total_records))
    n = 1;  %initialize num. reports in this event
    clear obs %working variable for each event
    %Initialize the output structure that we add to within each event and then
    %erase before the next event:
    obs.time(1) = total_records(index); %initialize
    obs.spd = zeros(1,length(stationnames));
    obs.spd(total_stationnums_of_records(index)) = total_spd(index); %"adds" windspeed to the current report's station. We'll divide it out by the station counts vector at the end of each event.
    obs.precip = zeros(1,length(stationnames));
    
    if total_precip(index) >= 0
        obs.precip(total_stationnums_of_records(index)) = total_precip(index); %assigns precip to the current report's station. ensures no -999 vals are included
    end
    obs.islightFZRA(1) = total_islightFZRA(index);
    obs.stationcounts = zeros(1,length(stationnames));
    obs.stationcounts(total_stationnums_of_records(index)) = 1;   %assign one count to whatever station this report is out.
    obs.stationlist(1) = total_stationnums_of_records(index);     %for nonevents
    
    %If next log is within max allowed nonevent hours, log it and keep
    %constructing an event. Nest in an if else to take care of the last
    %record.
    
    while(total_records(index + 1) <= ...
            total_records(index) + hours(max_nonevent_hrs))
        %Add this record's information to obs, the current event structure:
        obs.time(n+1) = total_records(index+1);
        obs.spd(total_stationnums_of_records(index+1)) = obs.spd(total_stationnums_of_records(index+1)) + total_spd(index+1); %"adds" windspeed to the current report's station. We'll divide it out by the station counts vector at the end of each event.
        if total_precip(index+1) >= 0
            obs.precip(total_stationnums_of_records(index+1)) = obs.precip(total_stationnums_of_records(index+1)) + total_precip(index+1); %assigns precip to the current report's station
        end
        obs.islightFZRA(n+1) = total_islightFZRA(index+1); %will just divide this at the end of the event by total num of reports.
        obs.stationcounts(total_stationnums_of_records(index+1)) = obs.stationcounts(total_stationnums_of_records(index+1)) + 1;   %assign one count to whatever station this report is from.
        obs.stationlist(n+1) = total_stationnums_of_records(index+1);    %Make a little list of the station numbers for use in nonevent counts.
        
        n = n+1;        %update number or reports within this event
        
        if index + 1 < length(total_records)    %(only look for more records if there are any remaining!)
            index = index + 1; %advance within total_records
        else %if we're at the last record, bail from this loop!
            break
        end
    end
        
    
    
    %Only include events with at least the min. number of reports.
    %If this event is included, record all its obs stuff:
    %Otherwise, skip them and don't add them to the output:
    if n >= min_reports_per_event
        %Round each observation to nearest timestep. This allows for use
        %directly with reanalysis or rawinsonde data:
        obs.time_unique = datenum(obs.time); %now in serial date, where unit is num of days.
        ts = timestep/24; %converted to serial date
        obs.time_unique = round(obs.time_unique/ts)*ts; %rounds to nearest timestep
        obs.time_unique = datetime(obs.time_unique,'ConvertFrom','datenum');    %convert back to datetime format
        obs.time_unique = unique(obs.time_unique); %eliminate duplicates
        
        %Now add it and the event ID to their respective vectors!
        event_times = [event_times obs.time_unique];
        event_ids = [event_ids id*ones(1,length((obs.time_unique)))];
        
        %Now calculate the event average variables!
        event_spd(id) = mean(obs.spd(obs.spd > 0)./obs.stationcounts(obs.spd > 0)); %so this is an average of the average windspeed at each station that participated in the event.
        event_spd_std(id) = std(obs.spd(obs.spd > 0)./obs.stationcounts(obs.spd > 0));
        %event_precip(id) = mean(obs.precip(obs.precip > 0));    %so this is an average of the total precip in hundredths of inches at each station that participated in the event.
        event_precip(id) = sum(obs.precip(obs.precip > 0).*(obs.stationcounts(obs.precip > 0)./sum(obs.stationcounts(obs.precip > 0)))); %this is a weighted average of total precip at a given station weighted by the stations that had the most reports
        event_precip_std(id) = std(obs.precip(obs.precip > 0).*(obs.stationcounts(obs.precip > 0)./sum(obs.stationcounts(obs.precip > 0))));
        event_pct_lightFZRA(id) = sum(obs.islightFZRA)/length(obs.islightFZRA);
        event_stationcounts(id,:) = obs.stationcounts;
        
        
        id = id + 1;    %increment event ID
        
    else
        %If this report doesn't meet the event threshold, add to the count of "nonevents" we get!
        nonevent_times = [nonevent_times obs.time];
        nonevent_stations = [nonevent_stations obs.stationlist];
        
    end
    
    %If there are any records left, move on to the next record, which will
    %be a new event.
    if index + 1 < length(total_records)    %(only look for more records if there are any remaining!)
        index = index + 1; %advance within total_records
    else %if we're at the last record, bail from this loop!
        break
    end
    
end


end

