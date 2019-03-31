%Script to download all daily surface pressure maps to an external
%harddrive.

clear

years = num2str([1979:2014]');

url_prmsl = [repmat('http://www.esrl.noaa.gov/psd/thredds/fileServer/Datasets/NARR/monolevel/prmsl.',36,1),years,repmat('.nc',36,1)];

filenames = url_prmsl(:,73:85);

for m = 1:length(filenames)
    tic
    outfilename = websave(filenames(m,:),url_prmsl(m,:))
    toc
    m
end
    here = websave('doggo','https://i.ytimg.com/vi/SfLV8hD7zX4/maxresdefault.jpg')