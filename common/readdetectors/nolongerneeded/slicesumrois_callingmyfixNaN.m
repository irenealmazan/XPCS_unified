function [sumorowsY,sumocolsX] = slicesumrois(imdata,ROIS,indxFLAG)
% [sum_over_rowsY,sum_over_columnsX] = slicesumrois(images,imageJ)
% Taking ROIs and summing values over 'short' and 'long' dimension to make image slices
%	Note - This does not 'normalize' the slices to the number of pixels summed
%			However, after using function, use sumorowY.image{i}./sumorowY.norm{i} 
%				for ROI {i} 
%
%			in calculations I refer to  Xcol or Yrow because
%				array ROWS are Y direction on image (imageJ or matlab reading)
%				array COLS are X direction on image (imageJ or matlab reading)
%					(there is no Xrow or Ycol)
%				NOTE THAT DURING BEAMLINE, SOMETIMES A ROTATION IS PUT ON
%					THE IMAGE FOR CONVENIENCE (BUT IT IS INCONVENIENT AFTERWARDS!)
%				IT DOES NOT CHANGE THE RAW DATA and note that the ROIS
%					that must be used in epics widget also refer to the RAW image
%
% CT 2016-10 - default uses Matlab indexing (ROIS given with 1,1 one corner)
%			   (this is old behavior)
%			   but if imageJ=1; it will use SPEC or IMAGEJ consistent indexing
%			   (ROI's with respect to 0,0 at one corner)
%					since imageJ is how ROI's are determined and set at beamline
% CT 2018-9  - note, sum(A,DIM,'omitnan') will sum over that dimension and 'omitnan's
%				replace with this command (and stop using my fixNaNimage)							
% 
% Output
%	[sumorowY,sumocolX]	 each are structures (each cell has different ROI)
%		within each structure, some things are saved as cells since they 
%		have different sizes per each ROI
%
%		sumorowY.image{1}, ...{2}.... will be the results of each ROI
%		sumorowY.ROIS  will be matrix and pass through matrix used for ROIs
%		sumorowY.ndx{1},..{2}  will be column vectors, and will be the index of unsummed region
%		sumorowY.indexing is string 'ImageJ ROI indexing assumed ...' or 'Matlab ROI indexing assumed'
%		sumorowY.norm  if want to normalize the images (image/(number of rows summed))
%						use sumorowY.image{1}./sumorowY.norm
%
%
% Input
%	images	(array)
%		2D or 3D matrix or structure with field {images} that is image matrix
%		3rd dimension is if there is a stack of images (with time or other)
%
%	ROIS
%		each row is an roi (lowX highX loY highY)
%		use showrois(ROIS) to show where the ROIS are on an image
%
%	indxFLAG (scalar)
%		= 1 (use imageJ type indexing, starting at [0:15] (16 pixels)
%		otherwise, use matlab type indexing, starting at [1:16] (16 pixels)
%
%
% CT 2018-03 work so that sum over X and Y works when those ROIS' have singleton 1 pixel
%		assume that 3rd dimension is never the singleton

% adminstrative homework
if nargin<2;	
	help slicesumrois;return;
end
if nargin<3;indxFLAG=0;end;

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

	ROIS = ROIS + indxFLAG;
	
if indxFLAG==1;
	% shifted spec indexing to work with matlab
	indexing = ['ImageJ ROI indexing assumed (0,0) start'];
else
	indexing = ['Matlab ROI indexing assumed (1,1) start'];
end

sumorowsY.indexing = indexing;
sumocolsX.indexing = indexing;

for ii=ndxrois

	colsndx = [ROIS(ii,1) : ROIS(ii,2)];
	rowsndx = [ROIS(ii,3) : ROIS(ii,4)];

	if indxFLAG==1	  % save as ImageJ but I already changed ROI lines
		sumorowsY.ndx{ii} = colsndx-1;
		sumocolsX.ndx{ii} = rowsndx-1;
	else
		sumorowsY.ndx{ii} = colsndx;
		sumocolsX.ndx{ii} = rowsndx;
	end

		sumorowsY.norm{ii}=length(colsndx);
		sumorowsX.norm{ii}=length(rowsndx);
% creates a column that it the length of the other index in ROI not summed over
% note, cannot transpose the XY if a 3D matrix, so must use permute 
% to effect the 2D transpose (sums are over columns)
%	sumorowsY.images{ii} = ...
%		squeeze( permute(sum(imdata(rowsndx,colsndx,:),'omitnan'),[2 1 3]) );
%	sumocolsX.images{ii} = ...
%		squeeze( permute(sum(permute(imdata(rowsndx,colsndx,:),[2 1 3]),'omitnan'),[2 1 3]) );

%size(imdata(rowsndx,colsndx,:))
% Fix for 1d lines CT 2018 03, replaced sum(fixNaN(...)  with sum(...,'omitnan') which is better
		[NY,NX,NZ]=size(imdata(rowsndx,colsndx,:));
		if all([NY NX]>1);

			sumorowsY.images{ii} = ...
				squeeze( permute(sum(            imdata(rowsndx,colsndx,:),'omitnan'),[2 1 3]) );
%				squeeze( permute(sum(fixNaNimage(imdata(rowsndx,colsndx,:))),[2 1 3]) );
			sumocolsX.images{ii} = ...
				squeeze( permute(sum(permute(            imdata(rowsndx,colsndx,:),[2 1 3]),'omitnan'),[2 1 3]) );
%				squeeze( permute(sum(permute(fixNaNimage(imdata(rowsndx,colsndx,:)),[2 1 3])),[2 1 3]) );
		
		else   % assume z is NOT singleton, but x or y is singleton but do need to treat differently
				% following should be done more elegantly, quick
			if NY==1  % already
				sumorowsY.images{ii} = squeeze(fixNaNimage(imdata(rowsndx,colsndx,:)));
				sumocolsX.images{ii} = squeeze(sum(fixNaNimage(imdata(rowsndx,colsndx,:),2)))';			
%				sumorowsY.images{ii} = squeeze(fixNaNimage(imdata(rowsndx,colsndx,:)));
%				sumocolsX.images{ii} = squeeze(sum(fixNaNimage(imdata(rowsndx,colsndx,:),2)))';
			
			elseif NX==1;
				sumorowsY.images{ii} = squeeze(sum(fixNaNimage(imdata(rowsndx,colsndx,:),1)))';
				sumocolsX.images{ii} = squeeze(fixNaNimage(imdata(rowsndx,colsndx,:)));	
%				sumorowsY.images{ii} = squeeze(sum(fixNaNimage(imdata(rowsndx,colsndx,:),1)))';
%				sumocolsX.images{ii} = squeeze(fixNaNimage(imdata(rowsndx,colsndx,:)));	
			end		
		end


end

end
%%%%%%%%%%%% HELPER FUNCTIONS AFTER THIS LINE %%%%%%%%%%%%
% function [Z]=helper1(inputts)
%
% end



