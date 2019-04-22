function [sumI,sumInorm] = sumrois(imdata,ROIS,indxFLAG)
% [Isum] = sumrois(images,ROIS,indxFLAG))
%		CT 2017-08 GBS has note that this doesn't work with single line ROIs  (fixed 2017-08CT
% Taking ROIs and summing values within, one number per ROI for each image in stack
% It handles correctly when images have NaN (i.e., it will ignore those points)
% 
% Output
%	Isum	Array  : each column for each ROI. Length of vectors (rows) num of images
%	sumInorm is array normed to the number of pixels in the ROI
% Input
%	images	2D or 3D matrix or structure with field {images} that is matrix
%	ROIS	each row is an roi (lowX highX loY highY)
%			use showrois(ROIS) to show where the ROIS are on an image
%			note for indexing pixels X is column and Y is row [Y,X]=size(ARRAY)
%   indxFLAG = 1; (=1 use ImageJ indexing) (= 0 [default] use matlab indexing)
%
% ROIS (assumed integer, as pixel) given as  Xmin Xmax Ymin Ymax on image
%
% CT 2016-10 - default uses Matlab indexing (ROIS given with 1,1 one corner)
%			   (this is old behavior)
%			   but if indxFLAG=1; it will use SPEC or IMAGEJ consistent indexing
%			   (ROI's with respect to 0,0 at one corner)
%					since imageJ is how ROI's are determined and set at 
%
% CT 2017-04 This has trouble if there are NaN in the plot(sums to NaN)
%				fixed by tweaking how it sums to avoid the NaN
%				these also are not counted in the normed sums
%
% CT 2017-08 Fixed so that can also use ROI's with only one pixel

if nargin<2;help sumrois;return;end
if nargin<3;indxFLAG = 0;end;

if isstruct(imdata);
	imdata = imdata.images;
end

ndxrois = [1:length(ROIS(:,1))];
%disp(['ndxrois ',int2str(ndxrois)]);

sizeROIS = [(ROIS(:,2)-ROIS(:,1)) .* (ROIS(:,4)-ROIS(:,3))];

if any([ROIS(:,2)<ROIS(:,1) | ROIS(:,4)<ROIS(:,3)])
	disp(['Expected ROIS are [Xstart Xend Ystart Yend], at least one of your ROIs has END smaller than START']);
	disp(['Behavior of program will be unexpected (it will just use the first point as a single point']);
end


for ii=ndxrois

		colsndx = [ROIS(ii,1) : ROIS(ii,2)]+indxFLAG;
		rowsndx = [ROIS(ii,3) : ROIS(ii,4)]+indxFLAG;

		% doing sums a bit more carefully to ignore the NaN that might be present in image
		ROIvecs		= reshape(imdata(rowsndx,colsndx,:),length(rowsndx).*length(colsndx),[]);
		ROIndx 		= ~isnan(ROIvecs);
			% replace NaN with zeros so sums work, ndx will be used for valid points
			ROIvecs(~ROIndx)=zeros(size(ROIvecs(~ROIndx)));
			
		sumI(:,ii) 		= sum(ROIvecs .* ROIndx,1)';
		sumInorm(:,ii)	= sumI(:,ii,1) ./ sum(ROIndx,1)';

%		sumI(:,ii) = squeeze(sum(sum(imdata(rowsndx,colsndx,:))));
%		sumInorm(:,ii) = sumI(:,ii)./sizeROIS(ii);

end

end
%%%%%%%%%%%% HELPER FUNCTIONS AFTER THIS LINE %%%%%%%%%%%%
