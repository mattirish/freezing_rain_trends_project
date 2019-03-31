%Plotting regional historical average frequencies with proxies for event frequencies.
regionalavgs = [16.0009850306071,9.73644271386903,20.8841057966377,15.5468655970494,10.3416082357219,14.8220033955296,21.1789767885304,19.3024325744809,14.7379051010770,15.0280760344692,13.0128955955620,8.07107115584570,7.81766687207378,15.8818772229840,18.9273669780661,17.6620086454449,10.0733966758025,17.6444048748030,16.4342682682604,19.9045263470005,12.7544587551426,15.4859301143409,13.4278061090542,14.7261226103160,10.9126453352780,13.7669953048176,14.1011502646740,14.1128103371281,11.1205342270259,16.2886740242043,10.4062942190775,23.1796064931857,23.1559381563711,13.9927670305075,7.92583946217283,14.8220033955296,14.8220033955296,19.9527676617467,10.0949016257389];
years = 1976:2014;

[P,S] = polyfit(years, regionalavgs,1);
[Y,delta] = polyval(P,years,S);


figure(1)
plot(years,regionalavgs)
hold on
plot(years,Y)
grid on
title('Regional Average Yearly Relative FZRA Frequency')
ylabel('Relative Frequency (% x 100) OR # Wiki Mentions')
xlabel('Year')
legend('Regional Avg.','39-yr Linear Best Fit')
hold on

%Add incidences of "freezing rain" from that year's Wiki article on
%"Global storm activity of [Year]". They're actually all US events.
load wikimentions
wikimentions(wikimentions==0) = NaN;
plot([1976:2014],wikimentions(:,2)')
legend('ASOS Rel. FZRA Freq.','Num. Wiki Article Mentions')