%FZRAmap_movie

%% 2D-Interpolating between points to make a raster map of average freezing rain freq
%We can use the griddata function to interpolate the data using
%the cubic method (better than linear interpolation for modeling
%things that vary smoothly across space)
%Initialize domain variables, etc.
resolution = 0.01;
x = (lon_lim(1):resolution:lon_lim(2));
y = (lat_lim(1):resolution:lat_lim(2));
[domain_y, domain_x] = meshgrid(x,y);

%Normalize to the total number of hours reported at each station:
numobservations = [numobsperyear_ia;
                    numobsperyear_il;
                    numobsperyear_in;
                    numobsperyear_mi;
                    numobsperyear_mn;
                    numobsperyear_ny;
                    numobsperyear_oh;
                    numobsperyear_onqb;
                    numobsperyear_pa;
                    numobsperyear_wi];

YearFreq = [a_ia.YearFreq;
                    a_il.YearFreq;
                    a_in.YearFreq;
                    a_mi.YearFreq;
                    a_mn.YearFreq;
                    a_ny.YearFreq;
                    a_oh.YearFreq;
                    a_onqb.YearFreq;
                    a_pa.YearFreq;
                    a_wi.YearFreq];
                
MonthFreq = [a_ia.MonthFreq;
                    a_il.MonthFreq;
                    a_in.MonthFreq;
                    a_mi.MonthFreq;
                    a_mn.MonthFreq;
                    a_ny.MonthFreq;
                    a_oh.MonthFreq;
                    a_onqb.MonthFreq;
                    a_pa.MonthFreq;
                    a_wi.MonthFreq];
              
%Find relative frequencies by dividing num FZRN hours recorded each year per station by number of hours recorded                
YearFreqrel = YearFreq./numobservations;
YearFreqrel = YearFreqrel*100*100; %scale by 100, then make percent

%Also include the four corner points of our region of interest with average
%frequencies so that the interpolation is carried out on the whole region:
stationlocations_wboundpts = [StationLocations; lat_lim' lon_lim'; lat_lim' fliplr(lon_lim)'];
MonthFreq_wboundpts = [MonthFreq; ones(1,12).*mean(MonthFreq); ones(1,12).*mean(MonthFreq); ...
    ones(1,12).*mean(MonthFreq); ones(1,12).*mean(MonthFreq)];
YearFreqrel_wboundpts = [YearFreqrel; ones(1,39).*mean(YearFreqrel); ones(1,39).*mean(YearFreqrel); ...
    ones(1,39).*mean(YearFreqrel); ones(1,39).*mean(YearFreqrel)];


%It's time! Interpolate the scattered data points using the griddata func:
%yearstoplot = 1976:2014; %select range of years to map average of
%yearstoplot = 1976:1990;
yearstoplot = 1990:2004;
%yearstoplot = 2005:2014;

%Set NaNs to the averages of their stations and figure out why they're NaNs
%later:
for k = 1:length(YearFreqrel_wboundpts)
    YearFreqrel_wboundpts(k,isnan(YearFreqrel_wboundpts(k,:))) = mean(YearFreqrel_wboundpts(k,~isnan(YearFreqrel_wboundpts(k,:))));
end



%------------------ BEGIN RECORDING LOOP -----------------------
vobj=VideoWriter('FZRNmovie1', 'MPEG-4');
vobj.FrameRate=2;
vobj.Quality=95
open(vobj);
for yearstoplot=1976:2014
    %Plot here:
    
    %Interpolate dat map:
    fzrnvals_wboundpts = mean(YearFreqrel_wboundpts(:,yearstoplot - 1975),2);
    fzrnmap = griddata(stationlocations_wboundpts(:,1),stationlocations_wboundpts(:,2), ...
        fzrnvals_wboundpts,domain_x,domain_y,'natural'); %beautiful! Last parameter describes the method of interpolation.
    
    
    %Plot the masked windmap! She's gonna be beautiful.
    figure(1)
    worldmap(lat_lim,lon_lim)
    colormap(parula)
    demcmap([min(fzrnmap),max(fzrnmap)],64,'window','window') %applies colormap to the raster data
    fzrn_ref = georasterref('RasterSize', size(fzrnmap), 'Latlim', lat_lim, 'Lonlim', lon_lim);
    geoshow(states,'FaceColor',[1 1 1])
    geoshow(provinces,'FaceColor',[1 1 1])
    %framem('ffacecolor',[.5 .7 .9]); %shows water as blue
    geoshow(fzrnmap,fzrn_ref,'DisplayType','texturemap','FaceAlpha',0.8)
    plotm(stationlocations_wboundpts(:,1),stationlocations_wboundpts(:,2),'r+') %plots measurement stations
    caxis([0, 25])
    colorbar
    hold on
    title('Relative Yearly Frequency of Hourly FZRN Reports')
    
    
    %Display the year. Add a timeline bar as a subplot later.
    mTextBox = uicontrol('style','text')
    set(mTextBox,'String',yearstoplot)
    % Something that I find useful is to set the Position Units to Characters, the default is pixels
    set(mTextBox,'Units','characters')
    % To move the the Text Box around you can set and get the position of Text Box itself
    mTextBoxPosition = get(mTextBox,'Position')
    % The array mTextBoxPosition has four elements
    % [x y length height]
    set(mTextBox,'Position',[85,.9,18,5])
    set(mTextBox,'FontSize',26)
    
    %Plot is finished. Get the frame and write it to the vid!
    F=getframe(gcf);
    writeVideo(vobj, F);
    cla(gca)
    yearstoplot
end
close(vobj)

