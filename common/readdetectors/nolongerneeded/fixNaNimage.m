function fixImage = fixNaNimage(Image,epilson)
% fixImage = fixNaNimage(Image)
%  NOTE  there is sum(...,NANFLAG) which does this better , 'omitnan'
% quick fix of NaN of image not ideal for handling NaN 
%   take NaN and turn to zero. 
%
% Ideally, in doing the sums for ROIS - I should eliminate these pixels
% and also eliminate them in calculating mean (which I do with
%  sumrois_wNaN which just sums over the total ROI)
% But I haven't worked on the slicesumrois yet (which is a bit more complicated)
%  ct notes 2017-04

if nargin<2;epsilon=0;end

fixImage = Image;  

%fixImage(find(Image<=0)) = epsilon.* ones(size(fixImage(find(Image<=0))));
%find(isnan(Image))

fixImage(find(isnan(Image))) = zeros(size(fixImage(find(isnan(Image)))));

