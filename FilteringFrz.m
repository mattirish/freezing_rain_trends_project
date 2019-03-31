clear

data = ['OriginalAlpena-726390         ';'OriginalBC-725396             ';'OriginalBentonH-726355        ';'OriginalDetroitColeman-725375 ';'OriginalDetroitWillow-725376  ';...
    'OriginalDTW-725370            ';'OriginalFlint-726370          ';'OriginalGR-726350             ';'OriginalHancock-727440        ';'OriginalHoughtonL-726380      ';'OriginalIronM-727437          ';...
    'OriginalJackson-725395        ';'OriginalKalamazoo-726357      ';'OriginalLansing-725390        ';'OriginalMuskegon-726360       ';'OriginalPellston-727347       ';...
    'OriginalPontiac-726375        ';'OriginalSaginaw-726379        ';'OriginalSaultSM-727340        ';'OriginalTraverseC-726387      '];
data = cellstr(data);


for i = 1:length(data)
    [index, double] = FrzRain(char(data(i)));
    totalindex(i, 1:length(index)) = index;
    totaldouble(i, 1:length(double)) = double;    
end
    
%%
clear
data = ['OriginalCarbondale_724336   ';'OriginalChampaign_725315    ';'OriginalChicagoDuPage_725305';'OriginalChicagoMidway_725340';...
    'OriginalChicagoOHare_725300 ';'OriginalDecatur_725316      ';'OriginalMoline_725440       ';'OriginalPeoria_725320       ';...
    'OriginalQuincy_724430       ';'OriginalRockford_725430     ';'OriginalSpringfield_724390  '];
data = cellstr(data);

for i = 1:length(data)
    [index, double] = FrzRain(char(data(i)));
    totalindex(i, 1:length(index)) = index;
    totaldouble(i, 1:length(double)) = double;    
end
    
   