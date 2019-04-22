function [sumI,sumInorm] = sumrois(imdata,ROIS,ImageJ)
% [Isum] = sumrois(images,ROIS,ImageJ))
% Taking ROIs and summing values within, one number per ROI for each image in stack
% 
% Output
%	Isum	Array  : each column for each ROI. Length of vectors (rows) num of images
%	sumInorm is array normed to the number of pixels in the ROI
% Input
%	images	2D or 3D matrix or structure with field {images} that is matrix
%	ROIS	each row is an roi (lowX highX loY highY)
%			use showrois(ROIS) to show where the ROIS are on an image
%			note for indexing pixels X is column and Y is row [Y,X]=size(ARRAY)
%   ImageJ = 1; (=1 assume ImageJ indexing in ROIS) = 0 (default, assume  matlab indexing used in ROIS
%
% ROIS (assumed integer, as pixel) given as  Xmin Xmax Ymin Ymax on image
%
% CT 2016-10 - default uses Matlab indexing (ROIS given with 1,1 one corner)
%			   (this is old behavior)
%			   but if indxFLAG=1; it will use SPEC or IMAGEJ consistent indexing
%			   (ROI's with respect to 0,0 at one corner)
%					since imageJ is how ROI's are determined and set at 
%
% CT 2017-08 need to Fix so that can also use ROI's with only one pixel

if nargin<2;help sumrois;return;end
if nargin<3;ImageJ = 0;end;


if isstruct(imdata);
	imdata = imdata.images;
end

ndxrois = [1:length(ROIS(:,1))];
%disp(['ndxrois ',int2str(ndxrois)]);

if any([ROIS(:,2)<ROIS(:,1) | ROIS(:,4)<ROIS(:,3)])

	disp(['Expected ROIS are [Xstart Xend Ystart Yend], at least one of your ROIs has END smaller than START']);
	disp(['Your ROIs with order error ',int2str(ROIS)]);
	disp(['Behavior of program will be unexpected (it will replace Xend or Yend to be same as start']);
	
		ND12  = find(ROIS(:,2)< ROIS(:,1));
		ND34  = find(ROIS(:,2) < ROIS(:,1));
			ROIS(ND12,2) = ROIS(ND12,1);
			ROIS(ND34,2) = ROIS(ND34,1);
	
	disp(['ROIs used to avoid order error (for now) in sumrois routine',int2str(ROIS)])
	
end


sizeROIS = [(ROIS(:,2)-ROIS(:,1)) .* (ROIS(:,4)-ROIS(:,3))];

for ii=ndxrois

		colsndx = [ROIS(ii,1) : ROIS(ii,2)]+ImageJ;
		rowsndx = [ROIS(ii,3) : ROIS(ii,4)]+ImageJ;

		
		% this appears to work for our regular ROIS, line ROIS, and point ROIS
		% In matlab R2006? the omitnan became available in sum and mean?
		sumI(:,ii) = squeeze(sum(sum(imdata(rowsndx,colsndx,:),1,'omitnan'),2,'omitnan'));
		sumInorm(:,ii) = squeeze(mean(mean(imdata(rowsndx,colsndx,:),1,'omitnan'),2,'omitnan'));


end

end
%%%%%%%%%%%% HELPER FUNCTIONS AFTER THIS LINE %%%%%%%%%%%%
