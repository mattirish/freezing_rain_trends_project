load geoid

figure(10000)
axesm eckert4;
framem; gridm;
axis off
geoshow(geoid, geoidrefvec, 'DisplayType', 'texturemap');
hcb = colorbar
set(get(hcb,'Xlabel'),'String','EGM96 Geoid Heights in Meters.')