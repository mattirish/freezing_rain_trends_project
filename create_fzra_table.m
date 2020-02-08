% Join metada for all FZRA observing stations to
% their respective stations in a.mat and b.mat

load a_all.mat %contains a, the batch summary data for all 98 stations
%load b_fixedcoords.mat %contains b, the full dataset of FZRA events for each station.
load b.mat %contains b, the full dataset of FZRA events for each station. *fixed b file. coords are fixed too. can prob delete b_fixedcoords

% Load trends data:
load yearly_trends.mat
load yearly_trends_duration.mat
load yearly_trends_intensity.mat

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

%%%%%%%%%%%% Output the new a and b from above if you're rerunning.

% Now to make that big summary table for the paper:
headers = {'State_Province','City','ICAO_ID','Median_Hours_per_Year','Mean_Hours_per_Year', ...
           'Decadal_Trend_Freq','p_value','Decadal_Trend_Freq_9614','p_value_9614', ...
           'Decadal_Trend_Dur','p_value_dur','Decadal_Trend_Dur_9614','p_value_dur_9614', ...
           'Decadal_Trend_Int','p_value_int','Decadal_Trend_Int_9614','p_value_int_9614'}
output_table2 = table(a.state, a.city, a.icao, median(a.YearFreq,2),mean(a.YearFreq,2), ...
                      10*yearly_trends_data(7,:)',yearly_trends_data(2,:)',10*yearly_trends_9614_data(7,:)',yearly_trends_9614_data(2,:)', ...
                      10*yearly_trends_duration_data(7,:)',yearly_trends_duration_data(2,:)',10*yearly_trends_duration_9614_data(7,:)',yearly_trends_duration_9614_data(2,:)', ...
                      10*yearly_trends_intensity_data(7,:)',yearly_trends_intensity_data(2,:)',10*yearly_trends_intensity_9614_data(7,:)',yearly_trends_intensity_9614_data(2,:)', ...
                      'VariableNames',headers)
writetable(output_table2,'summary_table_v2.csv')


% Create state average frequencies:
mean_freq_states = zeros(1,length(unique(a.state)));
median_freq_states = zeros(1,length(unique(a.state)));
all_states = unique(a.state);
for state = 1:length(all_states)
    mean_freq_states(state) = mean(a.YearFreq(strcmp(a.state,all_states(state))));
    median_freq_states(state) = median(a.YearFreq(strcmp(a.state,all_states(state))));
end

% Output state trends:
headers = {'State_Province','City','ICAO_ID','Median_Hours_per_Year','Mean_Hours_per_Year', ...
           'Decadal_Trend_Freq','p_value','Decadal_Trend_Freq_9614','p_value_9614', ...
           'Decadal_Trend_Dur','p_value_dur','Decadal_Trend_Dur_9614','p_value_dur_9614', ...
           'Decadal_Trend_Int','p_value_int','Decadal_Trend_Int_9614','p_value_int_9614'}
output_table3 = table(unique(a.state), repmat(cellstr(''),length(unique(a.state)),1), repmat(cellstr(''),length(unique(a.state)),1), median_freq_states', mean_freq_states', ...
                      10*str2double(yearly_trends_state_avg.('Sen''s slope')),str2double(yearly_trends_state_avg.('new P-value')),str2double(yearly_trends_state_avg_9614.('new P-value')),str2double(yearly_trends_state_avg_9614.('Sen''s slope')), ...
                      10*str2double(yearly_trends_duration_state_avg.('Sen''s slope')),str2double(yearly_trends_duration_state_avg.('new P-value')),10*str2double(yearly_trends_duration_state_avg_9614.('Sen''s slope')),str2double(yearly_trends_duration_state_avg_9614.('new P-value')), ...
                      10*str2double(yearly_trends_intensity_state_avg.('Sen''s slope')),str2double(yearly_trends_intensity_state_avg.('new P-value')),10*str2double(yearly_trends_intensity_state_avg_9614.('Sen''s slope')),str2double(yearly_trends_intensity_state_avg_9614.('new P-value')), ...
                      'VariableNames',headers)
writetable(output_table3,'summary_table_v2_states.csv')


% Output domain-wide trends:
headers = {'State_Province','City','ICAO_ID','Median_Hours_per_Year','Mean_Hours_per_Year', ...
           'Decadal_Trend_Freq','p_value','Decadal_Trend_Freq_9614','p_value_9614', ...
           'Decadal_Trend_Dur','p_value_dur','Decadal_Trend_Dur_9614','p_value_dur_9614', ...
           'Decadal_Trend_Int','p_value_int','Decadal_Trend_Int_9614','p_value_int_9614'}
output_table4 = table(cellstr('Domain Average'), cellstr(''), cellstr(''), median(median(a.YearFreq,2)),mean(mean(a.YearFreq,2)), ...
                      10*yearly_trends_domain_avg(7,:)',yearly_trends_domain_avg(2,:)',10*yearly_trends_domain_avg_9614(7,:)',yearly_trends_domain_avg_9614(2,:)', ...
                      10*yearly_trends_duration_domain_avg(7,:)',yearly_trends_duration_domain_avg(2,:)',10*yearly_trends_duration_domain_avg_9614(7,:)',yearly_trends_duration_domain_avg_9614(2,:)', ...
                      10*yearly_trends_intensity_domain_avg(7,:)',yearly_trends_intensity_domain_avg(2,:)',10*yearly_trends_intensity_domain_avg_9614(7,:)',yearly_trends_intensity_domain_avg_9614(2,:)', ...
                      'VariableNames',headers)
writetable(output_table4,'summary_table_v2_domain.csv')




