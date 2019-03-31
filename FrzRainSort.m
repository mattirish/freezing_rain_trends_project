function [] = FrzRainSort( FreezingData, long, double, num )
%FrzRainSort  Sorts previously filtered data by FrzRainCompile and outputs
% new .txt file
%   Needs the FreezingData matrix, length(index), double vector, and num indicator
%   from the previous function in that order as input

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creating a second txt file with the 
% data, 
% number of events/ hr,
% automated/manual

% Organizing Date
fid = fopen(FreezingData);
strdata = fgetl(fid);

% Reading the data from the txt file created in FrzRainCompile into matrices
for i = 1:long % long indicates the length of the txt file created in FrzRainCompile
    strdata = fgetl(fid);
    Year(i) = str2double(strdata(14:17));
    Month(i) = str2double(strdata(18:19));
    Day(i) = str2double(strdata(20:21));
    Hour(i) = str2double(strdata(22:23));
    YearMDH(i) = str2double(strdata(14:23)); % Full vector of year, month, day and hour
    FullDate(i) = str2double(strdata(14:25)); % Full vector of year, month, day, hour, minute 
    % Matrix contains the numbers listed under manual observations (rows
    % 1-4) and autmated observations (5-8), which simply indicate the type
    % of weather occuring at the corresponding time
    matrix(1,i) = str2double(strdata(58:59));
    matrix(2,i) = str2double(strdata(61:62));
    matrix(3,i) = str2double(strdata(64:65));
    matrix(4,i) = str2double(strdata(67:68));
    matrix(5,i) = str2double(strdata(70:71));
    matrix(6,i) = str2double(strdata(73:74));
    matrix(7,i) = str2double(strdata(76:77));
    matrix(8,i) = str2double(strdata(79:80));
end


% Finding number of occurences per hour 
EndLoop = FullDate(long);

i = 1;
z = 0;
while(FullDate(i) ~= EndLoop) 
    z = z + 1;
    index2(z) = i;
    Year2(z) = Year(i);
    Month2(z) = Month(i);
    Day2(z) = Day(i);
    Hour2(z) = Hour(i);
    count(z) = 1;
    k = 1;
    while(k ~= 0)
        if i >= length(YearMDH) % guarding against a possible attempt to pull from an out of bounds incex
            k = 0;
            i = long;
        elseif YearMDH(i) == YearMDH(i+k); % If the event has the same YearMDH, it happened in the same hour
            count(z) = k+1; % Vector of the number of observations that took place in the hour that defined each event
            k = k + 1; 
        else
            i = i + k;
            k = 0;
        end
        if (i+k) > length(YearMDH) % Again guarding against out of bounds indices
            i = i + (k-1);
            k = 0;
        end
    end
end

% Protecting against the possibility of a double observation occuring on
% the final line of the txt file from FrzRainCompile
if double(length(double)) == long
    count(length(count)) = count(length(count)) + 1;
end

% Finding method of measurement
% Getting rid of any number not associated with a freezing rain event
% (other types of observations can occur at the same time as a Freezing Rain observation)
% For Manual
for i = 1:4
    for j = 1:long
        if matrix(i,j) < 66 || matrix(i,j) > 67
            matrix(i,j) = 0;
        end
    end
end

% For Automated
for i = 5:8
    for j = 1:long
        if matrix(i,j) < 64 || matrix(i,j) > 66
            matrix(i,j) = 0;
        end
    end
end

for i = 1:length(count)
    auto = 0;
    manual = 0;
    for j = 1:count(i)
        z = index2(i) + (j-1); 
        if ((matrix(1,z) > 0) || (matrix(2,z) > 0) || (matrix(3,z) > 0) || (matrix(4,z) > 0))
            manual = 1;
        end
        if ((matrix(5,z)) > 0 || (matrix(6,z) > 0) || (matrix(7,z) > 0) || (matrix(8,z) > 0))
            auto = 1;
        end
    end
    if ((manual > 0) && (auto > 0))
        measure(i) = 3; % 3 if both auto and manual observation of freezing rain occurred within the hour
    elseif (manual > 0)
        measure(i) = 1; % 1 if only manual observations of freezing rain occurred within the hour
    elseif (auto > 0)
        measure(i) = 2; % 2 if only auto observations of freezing rain occurred within the hour
    end
end

% Writing filtered txt file
occurred(1,1:length(count)) = 'Y';     
str = sprintf('FreezingRainFilter_%i.txt',num);
fid = fopen(str, 'w');
fprintf(fid, '%4s %5s %3s %4s %5s %5s %4s\n', 'Year', 'Month', 'Day', 'Hour', '(Y/N)', 'Freq.', 'Type');
for i = 1:length(count)
    fprintf(fid,'%4i %3i %4i %4i %4s %5i %5i\n',Year2(i)',Month2(i)',Day2(i)',Hour2(i)',occurred(i)',count(i)',measure(i)');
end



end
