function fixImage = fixlogimage(Image,epilson)
% fixImage = fixlogimage(Image)
% take an Image, and all zeros and negative numbers make zero
% (the surf and pcolor can handle log(0) but not log(negatives)
%

if nargin<2;epsilon=1e-50;end

fixImage = Image;  

%fixImage(find(Image<=0)) = epsilon.* ones(size(fixImage(find(Image<=0))));

fixImage(find(Image<=0)) = zeros(size(fixImage(find(Image<=0))));

