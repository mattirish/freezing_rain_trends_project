clear
clc

%% Master Extraction file for all ASOS stations across states and provincesS
%Extract the fzra data needed into structures.

%Edit: 8 April 2018. Last 8 values of fields MW1 .. PCP01 are messed up.
%They're showing up as the value of the column to the left of the intended
%column in the filtered text file. Trying to fix this.
%Also adding AW1...AW4 fields for intensity considerations.

%%
%Paste station lats and lons into a stationlocations vector:
%All stations (updated from ASOS_Stations.xls on 8 Feb 2017):
%(Spaces are locations of thrown out stations in the alphabet)

%These station locations below are NOT right. The b-file as of 16 April has
%been fixed using the coords from b_fixedcoords.
% StationLocations = [46.49	-84.51
% 
% 42.267	-82.967
% 	
% 	
% 45.3225	-75.6692
% 44.75	-81.1
% 48.0564	-77.7867
% 48.2458	-79.0342
% 46.625	-80.7989
% 46.367	-79.417
% 45.8833	-82.5667
% 48.5697	-81.3767
% 48.367	-89.317
% 49.4139	-82.4675
% 40.19319067	-76.76262061
% 39.87224944	-75.24086583
% 	
% 39.99797222	-82.89188889
% 39.94444719	-81.89209403
% 	
% 39.90225139	-84.21941019
% 	
% 40.70747797	-84.02707806
% 38.03841667	-87.53086111
% 	
% 	
% 39.71729914	-86.29466119
% 40.41230556	-86.93688889
% 39.84420989	-89.67808981
% 39.94285269	-91.19458375
% 40.77725	-73.87261111
% 41.62658278	-73.8841925
% 41.06694444	-73.70755556
% 40.79525	-73.10022222
% 	
% 	
% 41.17827778	-78.89869444
% 40.29636111	-78.32002778
% 41.33847222	-75.72338889
% 42.20855556	-75.97972222
% 42.15988889	-76.89161111
% 40.65236267	-75.44040119
% 42.74911111	-73.80197222
% 43.11118694	-76.10631056
% 40.49147222	-80.23286111
% 40.91608333	-81.44219444
% 43.34122222	-73.61030556
% 41.40941669	-81.85498044
% 41.26073556	-80.67909667
% 42.08202139	-80.17621556
% 41.80306778	-78.64012083
% 42.94052472	-78.73216667
% 43.10733333	-78.94619444
% 43.11886111	-77.67238889
% 41.98164861	-87.90667142
% 41.90694386	-88.24820314
% 40.03883333	-88.27780556
% 39.8345625	-88.86568917
% 40.66420139	-89.69325778
% 40.97847222	-85.19516667
% 	
% 41.7833	-87.75
% 41.70822222	-86.31733333
% 41.58680556	-83.80783333
% 41.01202778	-83.66861111
% 42.30727806	-85.25147972
% 42.77863889	-84.58619444
% 42.26041667	-84.46044444
% 44.35980556	-84.67111111
% 42.66563619	-83.42050564
% 42.40919444	-83.00986111
% 40.783225	-91.12550556
% 42.19536111	-89.09722222
% 41.44828558	-90.50752011
% 41.88510158	-91.71228861
% 41.53397222	-93.66308333
% 41.10659611	-92.44793972
% 	
% 42.55708139	-92.40034361
% 43.1577925	-93.33126056
% 42.40260333	-96.38436694
% 44.93583333	-74.84555556
% 43.99192222	-76.02173861
% 	
% 42.12858333	-86.4285
% 	
% 43.53291667	-84.07963889
% 44.74163889	-85.58236111
% 42.88083333	-85.52280556
% 46.47922222	-84.36838889
% 45.57091667	-84.79672222
% 	
% 47.16841722	-88.48906083
% 42.21244444	-83.35338889
% 42.94721072	-87.89672956
% 43.13985778	-89.33751361
% 	
% 43.87926556	-91.25663803
% 44.86580453	-91.48425517
% 43.90827778	-92.50002778
% 44.48463889	-88.12972222
% 	
% 44.92628453	-89.62700175
% 	
% 	
% 	
% 	
% 44.54688889	-95.082
% 45.86630556	-95.39466667
% 44.88195667	-93.22176556
% 45.81836111	-88.11455556
% 45.07808333	-83.56030556
% 	
% 	
% 42.23438889	-85.55155556
% 46.84208333	-92.19363889
% 47.38658333	-92.83897222
% 48.56558333	-93.40216667
% 40.63975111	-73.77892556];
StationLocations = [46.49	-84.51
42.267	-82.967
45.3225	-75.6692
44.75	-81.1
48.0564	-77.7867
48.2458	-79.0342
46.625	-80.7989
46.367	-79.417
45.8833	-82.5667
48.5697	-81.3767
48.367	-89.317
49.4139	-82.4675
40.19319067	-76.76262061
39.87224944	-75.24086583
39.99797222	-82.89188889
39.94444719	-81.89209403
39.90225139	-84.21941019
40.70747797	-84.02707806
38.03841667	-87.53086111
39.71729914	-86.29466119
40.41230556	-86.93688889
39.84420989	-89.67808981
39.94285269	-91.19458375
40.77725	-73.87261111
41.62658278	-73.8841925
41.06694444	-73.70755556
40.79525	-73.10022222
41.17827778	-78.89869444
40.29636111	-78.32002778
41.33847222	-75.72338889
42.20855556	-75.97972222
42.15988889	-76.89161111
40.65236267	-75.44040119
42.74911111	-73.80197222
43.11118694	-76.10631056
40.49147222	-80.23286111
40.91608333	-81.44219444
43.34122222	-73.61030556
41.40941669	-81.85498044
41.26073556	-80.67909667
42.08202139	-80.17621556
41.80306778	-78.64012083
42.94052472	-78.73216667
43.10733333	-78.94619444
43.11886111	-77.67238889
41.98164861	-87.90667142
41.90694386	-88.24820314
40.03883333	-88.27780556
39.8345625	-88.86568917
40.66420139	-89.69325778
40.97847222	-85.19516667
41.7833	-87.75
41.70822222	-86.31733333
41.58680556	-83.80783333
41.01202778	-83.66861111
42.21244444	-83.35338889
42.40919444	-83.00986111
42.23930367	-83.53096528
42.77863889	-84.58619444
42.26041667	-84.46044444
42.30727806	-85.25147972
40.783225	-91.12550556
42.19536111	-89.09722222
41.44828558	-90.50752011
41.88510158	-91.71228861
41.53397222	-93.66308333
41.10659611	-92.44793972
42.55708139	-92.40034361
43.1577925	-93.33126056
42.40260333	-96.38436694
44.93583333	-74.84555556
43.99192222	-76.02173861
42.88083333	-85.52280556
42.23438889	-85.55155556
43.16767172	-86.23543872
42.96541667	-83.74366667
42.66563619	-83.42050564
43.53291667	-84.07963889
44.74163889	-85.58236111
45.07808333	-83.56030556
42.94721072	-87.89672956
43.13985778	-89.33751361
43.87926556	-91.25663803
44.86580453	-91.48425517
43.90827778	-92.50002778
44.48463889	-88.12972222
44.92628453	-89.62700175
44.54688889	-95.082
45.86630556	-95.39466667
44.88195667	-93.22176556
46.47922222	-84.36838889
45.57091667	-84.79672222
47.16841722	-88.48906083
46.84208333	-92.19363889
47.38658333	-92.83897222
48.56558333	-93.40216667
40.63975111	-73.77892556];


%% Make a structure of all the filenames in this directory that end in ".txt".
%They should be in the alphabetical order to match StationLocations above.
filenames = dir('*.txt');

YearFreq = zeros(length(filenames),39);
MonthFreq = zeros(length(filenames),12);

%We'll have two structures:  1) a summary structure with matrices for all
                            %stations and data, (a_allstations_batch) and
                            %2) a struct containing a struct with full data
                            %for each station (b_all)
                            
count = 0;                           
for q = 1:length(filenames)
    count = count+1;
    disp(count)
    
    fid = fopen(filenames(q).name,'r');
    %fid = fopen('FZRN_725370.txt','r');
    %Format reference: ftp://ftp.ncdc.noaa.gov/pub/data/noaa/ISH-DVD2012/ish-abbreviated.txt
    %Format:
    %      USAF YR  MO  DA  HR  DIR  SPD  MW1  MW2  MW3  MW4  TEMP DEWP  SLP  PCP01
    % Data:   1  2   3   4   5    6    7    8    9   10   11    12   13   14     15
    %      USAF YR  MO  DA  HR  DIR  SPD  MW1  MW2  MW3  MW4  AW1  AW2  AW3  AW4  TEMP  DEWP  SLP  PCP01
    % Data:   1  2   3   4   5    6    7    8    9   10   11   12   13   14   15    16    17   18     19                                   
    dataformat = '%f %*f %f %f %f %f %*f %f %f %*f %*f %*4c %*f %*f %*f %*f %f %f %f %f %f %f %f %f %*f %f %f %f %*f %*f %*f %*f %f %*f %*f %*f %*f %*[^\n]';
    Data = textscan(fid, dataformat,'HeaderLines',1,'Delimiter' , '\r','MultipleDelimsAsOne',0); %added that last 'MultipleDelimsAsOne',0 to fix the issue. It was reading OVC as 3 characters and then got messed up when it got to -999.
    
    
    stn = cell2mat(Data(1));
    stn = num2str(stn(1));
    stn = strcat('USAF',stn);
    b_all.(stn).YR = cell2mat(Data(2));
    b_all.(stn).MO = cell2mat(Data(3));
    b_all.(stn).DA = cell2mat(Data(4));
    b_all.(stn).HR = cell2mat(Data(5));
    b_all.(stn).DIR = cell2mat(Data(6));
    b_all.(stn).SPD = cell2mat(Data(7));
    b_all.(stn).MW1 = cell2mat(Data(8));
    b_all.(stn).MW2 = cell2mat(Data(9));
    b_all.(stn).MW3 = cell2mat(Data(10));
    b_all.(stn).MW4 = cell2mat(Data(11));
    b_all.(stn).AW1 = cell2mat(Data(12));
    b_all.(stn).AW2 = cell2mat(Data(13));
    b_all.(stn).AW3 = cell2mat(Data(14));
    b_all.(stn).AW4 = cell2mat(Data(15));
    b_all.(stn).TEMP = cell2mat(Data(16));
    b_all.(stn).DEWP = cell2mat(Data(17));
    b_all.(stn).SLP = cell2mat(Data(18));
    b_all.(stn).PCP01 = cell2mat(Data(19));
    
    %Other key fields for each station:
    b_all.(stn).coords = StationLocations(q,:);
    
    
    %Now for the summary structure (still making it for legacy purposes):
    %Record average wind speed and direction for each station:
    Dir = b_all.(stn).DIR;
    WindDir(q) = median(Dir(Dir >= 0));
    Speed = b_all.(stn).SPD;
    WindSpeed(q) = mean(Speed(Speed >= 0));
    
    % Monthly Proportions? Monthly Frequencies (Month Occurence/(total hours in Month))
    % Month(1) corresponds to January, etc
    for j = 1:12
        % Finding total freezing rain events over full time span. (This has
        % to be normalized for total num. observations, so it's inaccurate.
        % Don't use this in final product--use b struct instead.
        MonthFreq(q,j) = sum(b_all.(stn).MO == j); % Rows is Stations, Columns is Months (Jan-Dec)
    end
    
    % Yearly Frequencies
    % Year(1) corresponds to 1976, Year(2) to 1977... Year(39) to 2014
    for k = 1:39
        YearFreq(q,k) = sum(b_all.(stn).YR == (1975 + k)); % Rows is Stations, Columns is Years (1976-2014)
        
        %Find hourly frequencies as a time series by nesting a loop here:
        for m = 1:24
            HourlyFreq_yearly(q,k,m) = sum(b_all.(stn).HR(b_all.(stn).YR == (1975 + k)) == m - 1); % dims are station, year, hour
        end
        
        %Find monthly frequencies as a time series by nesting a loop here:
        for m = 1:12
            MonthFreq_yearly(q,k,m) = sum(b_all.(stn).MO(b_all.(stn).YR == (1975 + k)) == m); % dims are station, year, month
        end
    end
    
    %Record hourly frequencies for hourly freq chart.
    for k = 1:24
        HourlyFreq(q,k) = sum(b_all.(stn).HR == (k - 1)); % Rows is Stations, Columns is Hours (1976-2014)
    end
    
    clear Format fid Data
end




%%
%Save this as a struct specific to the state for input in
%Matt_spatialtrends.m:
a_allstations_batch = struct('YearFreq',YearFreq, ...
            'HourlyFreq',HourlyFreq, ...
            'HourlyFreq_yearly',HourlyFreq_yearly, ...
            'MonthFreq',MonthFreq, ...
            'MonthFreq_yearly',MonthFreq_yearly, ...
            'WindSpeed',WindSpeed, ...
            'WindDir',WindDir, ...
            'StationLocations',StationLocations);

save 'a_allstations_batch.mat' a_allstations_batch
save 'b_all.mat' b_all

