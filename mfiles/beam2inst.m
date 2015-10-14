function y=beam2inst(ba,c)
% Creates transformation matrix for converting from beam coordinates to 
% instrument coordinates based the beam angle.  
%
% ba is the beam angle of the ADCP (typically either 20 or 30 degrees)
% c is whether the instrument is convex (c=1) or concave (c other than 1). Nearly
% all ADCPs are convex (some early vessel-mounted systems were concavbe).
%
% Create the apropriate scale factors:
%
barad=ba*pi/180;
a=1/(2*sin(barad));
b=1/(4*cos(barad));
d=a/sqrt(2);
%
% Create the matrix:
%
y=[c*a -c*a 0 0
   0 0 -c*a c*a
   b b b b
   d d -d -d];
return