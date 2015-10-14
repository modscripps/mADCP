function y=inst2earth(head,pitch,roll,updown)
% Creates transformation matrix for converting from instrument coordinates to 
% earth coordinates based on heading, pitch and roll information.  Assumes
% arguments are given in degrees.  Updown should be 1 for a down-looking system,
% and some number other than 1 for an up-looking system.
%
% Note that the matrix created is 4x4, with the additional row and column added
% so that the error velocity calculated in beam2inst.m will carry through.
%
% Convert input degrees to radians
hrad=head*pi/180;
prad=pitch*pi/180;
if updown == 1
   rrad=roll*pi/180;
else
   rrad=(roll+180)*pi/180;
end
%
% Create the heading rotation matrix:
%
h=[cos(hrad)  sin(hrad)  0 0
   -sin(hrad) cos(hrad) 0 0
   0 0 1 0
	0 0 0 1];
%
% Create the pitch rotation matrix:
%
p=[1 0 0 0
   0  cos(prad) -sin(prad) 0
   0  sin(prad) cos(prad) 0
	0 0 0 1];
%
% Create the roll rotation matrix:
%
r=[cos(rrad) 0 sin(rrad) 0
   0  1  0 0
   -sin(rrad) 0 cos(rrad) 0
	0 0 0 1];
%
% create the final matrix
%
y=h*p*r;
return