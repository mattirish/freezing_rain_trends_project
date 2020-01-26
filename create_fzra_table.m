% Join metada for all FZRA observing stations to
% their respective stations in a.mat and b.mat

load a_all.mat %contains a, the batch summary data for all 98 stations
%load b_fixedcoords.mat %contains b, the full dataset of FZRA events for each station.
load b.mat %contains b, the full dataset of FZRA events for each station. *fixed b file. coords are fixed too. can prob delete b_fixedcoords


metadata = readtable('observing_stations_metadata.csv');


%Make assignments to a. (Everything's in the right order already)
%(It's bad practice to rely on the order of a struct but it's fine):
test = table2cell(metadata);
a.('icao')      = metadata.x___icao;
a.('city')      = metadata.city;
a.('state')     = metadata.state;
a.('id')        = metadata.id;
a.('lat')       = metadata.lat;
a.('lon')       = metadata.lon;

names = fieldnames(b);
for station = 1:length(fieldnames(b))
    b.(names{station}).('icao')      = metadata.x___icao(station);
    b.(names{station}).('city')      = metadata.city(station);
    b.(names{station}).('state')     = metadata.state(station);
    b.(names{station}).('id')        = metadata.id(station);
    if b.(names{station}).coords == [metadata.lat(station),metadata.lon(station)]
        'uh oh'
        break
    end
    b.(names{station}).('lat')       = metadata.lat(station);
    b.(names{station}).('lon')       = metadata.lon(station);
    
end

%%%%%%%%%%%%Output the new a and b from above if you're rerunning.

% Now to make that big summary table for the paper:
headers = {'State_Province','City','ICAO_ID','Median_Hours_per_Year','Mean_Hours_per_Year','Decadal_Trend','p_value'}
output_table2 = table(a.state, a.city, a.icao, median(a.YearFreq,2),mean(a.YearFreq,2),10*yearly_trends_data(7,:)',yearly_trends_data(2,:)','VariableNames',headers)
writetable(output_table2,'summary_table_2010.csv')
