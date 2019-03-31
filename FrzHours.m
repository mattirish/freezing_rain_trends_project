function [ hours ] = FrzHours( DATA )
%FrzHours Finds how many hours of data available in a certain type of data set
%   

Format = {'%13c%12c%121c','HeaderLines',1,'Delimiter' , '\r\n'};
fid = fopen(DATA,'r');
Data = textscan(fid, Format{:});
date = cell2mat(Data(2));
hours = 1;

for i = 1:(length(date)-1)
    if str2double(date(i,9:10)) ~= str2double(date(i+1,9:10))
        hours = hours + 1;
    end
end



end

