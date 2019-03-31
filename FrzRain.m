function [index, double] = FrzRain( FreezingData )
%FrzRain Combines two fucntions, FrzRainCompile and FrzRainSort 
%   Outputs two txt files with the desired freezing rain data 
%   based on a unique input of data

[index, double, num, str] = FrzRainCompile(FreezingData);

long = length(index);

FrzRainSort(str,long,double,num);


end

