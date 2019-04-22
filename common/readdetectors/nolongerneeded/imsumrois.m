function [sumorowsY,sumocolsX] = slicesumrois(imdata,ROIS)
% [sum_over_rowsY,sum_over_columnsX] = slicesumrois(images,ROIS)
%		I think this one isn't used CT 2018-08
% Taking ROIs and summing values over 'short' and 'long' dimension to make image slices
%	Note - typically we are interested in sum over short length in ROI so 
%			if one really wants summed over long, then ask for it in output
%			in calculations I refer to  Xcol or Yrow because
%				array ROWS are Y direction on image
%				array COLS are X direction on image
%					(there is no Xrow or Ycol)
%							
% 
% Output
%	[sumorowY,sumocolX]	 each are structures (each cell has different ROI)
%		within each structure, some things are saved as cells since they 
%		have different sizes per each ROI
%
%		sumorowY.image{1}, ...{2}.... will be the results of each ROI
%		sumorowY.ROIS  will be matrix and pass through matrix used for ROIs
%		sumorowY.ndx{1},..{2}  will be column vectors, and will be the index of unsummed region

%
% Input
%	images	(array)
%		2D or 3D matrix or structure with field {images} that is image matrix
%		3rd dimension is if there is a stack of images (with time or other)
%
%	ROIS
%		each row is an roi (lowX highX loY highY)
%		use showrois(ROIS) to show where the ROIS are on an image

% adminstrative homework
if nargin<2;	
	help sumrois;return;
end
if isstruct(imdata);
	imdata = imdata.images 
end

Ncell = length(ROIS(:,1));

sumorowsY.images = cell(1,Ncell);
sumocolsX.images = cell(1,Ncell);
sumorowsY.ndx = cell(1,Ncell);
sumocolsX.ndx = cell(1,Ncell);


ndxrois = [1:Ncell];

	sumorowsY.ROIS = ROIS;
	sumocolsX.ROIS = ROIS;

for ii=ndxrois

	colsndx = [ROIS(ii,1) : ROIS(ii,2)];
	rowsndx = [ROIS(ii,3) : ROIS(ii,4)];
	
	sumorowsY.ndx{ii} = colsndx;
	sumocolsX.ndx{ii} = rowsndx;

% creates a column that it the length of the other index in ROI not summed over
% note, cannot transpose the XY if a 3D matrix, so must use permute 
% to effect the 2D transpose (sums are over columns)
	sumorowsY.images{ii} = ...
		squeeze( permute(sum(imdata(rowsndx,colsndx,:)),[2 1 3]) );
	sumocolsX.images{ii} = ...
		squeeze( permute(sum(permute(imdata(rowsndx,colsndx,:),[2 1 3])),[2 1 3]) );

end

end
%%%%%%%%%%%% HELPER FUNCTIONS AFTER THIS LINE %%%%%%%%%%%%
% function [Z]=helper1(inputts)
%
% end



