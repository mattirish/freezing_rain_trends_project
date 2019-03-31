%Generate's movie of a randomly moving circle
%From http://stackoverflow.com/questions/14949507/how-to-generation-animation-as-mpg-or-speed-up-gif-in-matlab


vobj=VideoWriter('MyMovieFile2', 'MPEG-4');
vobj.FrameRate=1;
vobj.Quality=75
open(vobj);
for i=1:100
  plot(rand,rand,'o')
  F=getframe(gcf);
  writeVideo(vobj, F);
  cla(gca)
end
close(vobj)