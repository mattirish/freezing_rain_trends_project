function [ available,frequency_year,frequency_hour_ava,frequency_hour_tot ] = FrzFreq( DATA, hours )
%FrzFreq Used in conjunction with FrzHours to find the frequency of freezing rain events as well as finds the consistency of measurement of the data set being used (certain type of data set)

Format = {'%35c','HeaderLines',1,'Delimiter' , '\r\n'};
fid = fopen(DATA,'r');
Data = textscan(fid, Format{:});
stuff = cell2mat(Data(1));
incidents = length(stuff); % Simply the length of the data set found in the Filtered Data (from FrzRainSort) txt file minus the header

frequency_year = incidents/39; % freezing rain events per year
frequency_hour_ava = incidents/hours; % freezing rain events out of the total available hours; freezing rain events per hour
frequency_hour_tot = incidents/341880; % freezing rain events out of the total possible hours

available = hours/341880; % available hours out of total possible hours



end

