function [index,double,num,str] = FrzRainCompile( DATA )
%FrzRainCompile Sorts unique format and creates a txt file with all freezing rain events
%   Outputs: index vector, index of doubles,
%   station number, and txt title string



%% Taking each  recorded freezing rain event and writing them into a txt file
fid = fopen(DATA);
strdata = fgetl(fid); % Getting rid of header line
i = 1;
k = 0;
num = str2num(DATA((length(DATA)-5):length(DATA)));
FreezingData = zeros(1500,147);
index = 0;
while(strdata ~= -1)
    i = i + 1;
    strdata = fgetl(fid);
    if strdata == -1
        strdata = -1;
    elseif str2double(strdata(58:59)) == 66 || str2double(strdata(58:59)) == 67
        k = k + 1;
        L = 1:length(strdata);
        index(k) = i; % Index of each line that has a freezing rain event
        FreezingData(k,L) = strdata; % Matrix created with each line of data with a freezing rain event
    elseif str2double(strdata(61:62)) == 66 || str2double(strdata(61:62)) == 67
        k = k + 1;
        L = 1:length(strdata);
        index(k) = i;
        FreezingData(k,L) = strdata;
    elseif str2double(strdata(64:65)) == 66 || str2double(strdata(64:65)) == 67
        k = k + 1;
        L = 1:length(strdata);
        index(k) = i;
        FreezingData(k,L) = strdata;
    elseif str2double(strdata(67:68)) == 66 || str2double(strdata(67:68)) == 67
        k = k + 1;
        L = 1:length(strdata);
        index(k) = i;
        FreezingData(k,L) = strdata;
    elseif str2double(strdata(70:71)) == 64 || str2double(strdata(70:71)) == 65 || str2double(strdata(70:71))== 66
        k = k + 1;
        L = 1:length(strdata);
        index(k) = i;
        FreezingData(k,L) = strdata;
    elseif str2double(strdata(73:74)) == 64 || str2double(strdata(73:74)) == 65 || str2double(strdata(73:74))== 66
        k = k + 1;
        L = 1:length(strdata);
        index(k) = i;
        FreezingData(k,L) = strdata;
    elseif str2double(strdata(76:77)) == 64 || str2double(strdata(76:77)) == 65 || str2double(strdata(76:77))== 66
        k = k + 1;
        L = 1:length(strdata);
        index(k) = i;
        FreezingData(k,L) = strdata;
    elseif str2double(strdata(79:80)) == 64 || str2double(strdata(79:80)) == 65 || str2double(strdata(79:80))== 66
        k = k + 1;
        L = 1:length(strdata);
        index(k) = i;
        FreezingData(k,L) = strdata;
    end
end

%check = FreezingData; % check will be output so it can be checked with the outputted txt file

FreezingData2 = FreezingData(1:length(index),1:147);
clear FreezingData

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Getting rid of any duplicates

k =1;
double(k) = 0;
for i = 2:length(index)
    if str2double(FreezingData2(i,14:25)) == str2double(FreezingData2((i-1),14:25))
        k = k + 1;
        double(k) = i+1; % index of each duplicate, will be output
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Outputting the first filtered freezing rain events txt file
str = sprintf('FreezingRain%i.txt',num);
fclose(fid);
fid = fopen(str, 'w');
fprintf(fid, '  USAF  WBAN YR--MODAHRMN DIR SPD GUS CLG SKC L M H  VSB MW MW MW MW AW AW AW AW W TEMP DEWP    SLP   ALT    STP MAX MIN PCP01 PCP06 PCP24 PCPXX SD\n');
for i = 1:length(index)
    fprintf(fid,'%s\n',FreezingData2(i,1:147)');
end







end

