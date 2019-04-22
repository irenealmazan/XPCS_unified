ANGLES = get_geoAngles(outdata,collabels,scninfo)
% ANGLES = get_geoAngles(outdata,collabels,scninfo)
% Typically for detector or hkl calculations, need all the geo angles
% Some may be in scan, and some may be in header (unmoved)
%  scninfo.geoangles_i are the angles at start of scan
%  scninfo.geoangles_labels are the names
%  outdata

% finds which columns in scan have motors that are geo motors
[NDX,NDXgeo] = chan2col(collabels,scninfo.geoangles_label,'quiet');


% but don't actually know which motors those are


ANGLESscn = 
