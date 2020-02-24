% Let's process and pull maps for all combinations of FZRA event definitions as a sensitivity analysis.
%Keep the system from sleeping while this runs:
system('caffeinate -dims &');

for min_reports_per_event = [3 5 6 10] %4:-1:2
    for max_nonevent_hrs = 12 %3:-1:1
        
        timestep = 3;               %rounds to every three hours
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
        clear dates
        
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
        %% Now apply the minimum per-station reports threshold to define the number
        %of reports needed to classify it as an ice storm and then
        %create a companion vector that classifies each event by index in increasing order.
        
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
            
            
            
            %Only include events with at least the min. number of reports at at least one station.
            %If this event is included, record all its obs stuff:
            %Otherwise, skip them and don't add them to the output:
            if any(obs.stationcounts >= min_reports_per_event)
            %if n >= min_reports_per_event
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
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Save the event times:
        A = exist('nonevents');
        if A == 0
            nonevents = 0;
        end
        save(strcat('eventtimes_v3_',num2str(min_reports_per_event),num2str(max_nonevent_hrs)), ...
            'event_ids','event_pct_lightFZRA','event_precip','event_precip_std',...
            'event_spd','event_spd_std','event_stationcounts','event_times','nonevents');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Pull down the maps:
        load latslons_subregion_for_PCA
        load monthlies_correctedfordomainsize %archive of MSL maps for all 456 months and 12 averages, one for each month!
        
        clear prmsl_anom
        clear airtemp
        clear hgt850_anom
        clear hgt1000_anom
        clear hgt500_anom
        clear prmsl_anom_mb
        clear hgt1000500_anom
                  
        % For dynamics analysis with NARR, we only want events from 1979
        % onward, so exclude the 1976-1978 events:
        event_times_case = event_times(event_times >= datetime(1979,1,1) & event_times < datetime(2015,1,1));
        event_ids = event_ids(event_times >= datetime(1979,1,1) & event_times < datetime(2015,1,1));
        
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
            n = 0;  %counting vector for num. 3-hrly reports in each event.
            
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
                while(event_ids(m+1) == event_ids(m) && m+1 < length(event_ids)) %While the next record is still in this event:
                    m = m + 1;
                    n = n + 1;
                end
            end
            
%             %Now find the index in the prmsl files to download.
%             %We'll use the map closest to the median of the event times.
%             dates(iter) = event_times_case(floor(m-(n)/2)); %track times corresponding with maps. This is the median-time of the event.
             
            %Or download the map closest to the hour in which the most FZRA
            %observations occur:
            dates(iter) = mode(event_times_case(event_ids == event_ids(m)));
            
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
            hgt850_anom(:,:,iter) = nc_varget(url_hgt,'hgt',[hgt_indextoplot 6 yinds xindw],[1 1 yindn xinde]); %format is (time,level,y,x)
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        display('printing maps. Value of first pixel on this iteration:')
        prmsl_anom(1,1,1)
        % Save the maps:
        save(strcat('maps_v3_',num2str(min_reports_per_event),num2str(max_nonevent_hrs)), ...
            'airtemp','dates','hgt1000500_anom','prmsl_anom_mb','hgt850_anom',...
            'event_times_case');
        
        %clear
    end
end

