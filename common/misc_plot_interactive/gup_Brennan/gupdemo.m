% gupdemo.m
% script demonstrating the use of the graphics tool 
% gup for both 2D and 3D plot
%
% The gup-package consist of:
% GUP.M CURSMODE.M GUPZOOM.M GUPINIT.M GUPINIT3.M GUPLMTS.M INBOUNDS.m
%
% Note:
% The global GUP_BUFF is an array of the datapoint
% clicked by the cursor for 2D plots.

% Copyright 1996 Anneli Munkholm & Sean M. Brennan.
% Stanford Synchrotron Radiation Laboratory
% Stanford Linear Accelerator Center, Stanford CA 94309
% Bren@slac.stanford.edu


figure
load logo
surf(L,R)
colormap(M)
gup
figure
x=[-1:.01:2];
plot(x+1,humps(x)+10)
gup
