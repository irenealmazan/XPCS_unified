function [sumorowsY,sumocolsX] = slicesumrois(imdata,ROIS,ImageJ)
% [sum_over_rowsY,sum_over_columnsX] = slicesumrois(images,ImageJ)
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
%				NOTE - ROIS are assumed to be given as 'indices' of araay 
%						not in other units e.g., angle or q
% CT 2018-9  - note, sum(A,DIM,'omitnan') will sum over that dimension and 'omitnan's
%				replace with this command (and stop using my fixNaNimage)
%				this is option in more recent versions of Matlab (newer than 2010?)	
%				also - in Matlab more recent than 2015 or 2016 log and log10 will handle NaN and 0
%				and pcolor and image and imagesc will handle arrays with log10(NaN or 0)
%				so I've taken out my fixes for if image has NaN or zero in it				
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
%	ImageJ (scalar)
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
if nargin<3;ImageJ=0;end;

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

	ROIS = ROIS + ImageJ;
	
if ImageJ==1;
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

	if ImageJ==1	  % save as ImageJ but I already changed ROI lines
		sumorowsY.ndx{ii} = colsndx-ImageJ;
		sumocolsX.ndx{ii} = rowsndx-ImageJ;
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
	%		'omitnan' flag wasn't available in matlab when some of these functions first were being created
		[NY,NX,NZ]=size(imdata(rowsndx,colsndx,:));
		if all([NY NX]>1);   % the usual type ROI's

			sumorowsY.images{ii} = ...
				squeeze(permute(sum(imdata(rowsndx,colsndx,:),1,'omitnan'),[2 1 3]) );
%				squeeze( permute(sum(fixNaNimage(imdata(rowsndx,colsndx,:))),[2 1 3]) );
			sumocolsX.images{ii} = ...
				squeeze(permute(sum(imdata(rowsndx,colsndx,:),2,'omitnan'),[2 1 3]) );
%				squeeze( permute(sum(permute(fixNaNimage(imdata(rowsndx,colsndx,:)),[2 1 3])),[2 1 3]) );
		
		elseif   any([NY NX]>1)    % if a line but not a point
				% while we assume z is NOT singleton,  x or y could be singleton but do need to treat differently
				% following should be done more elegantly, quick, 
				%  if single point ROI this will not work!! (still need to fix in case someone asks to 'sum'
			if NY==1  % already
				sumorowsY.images{ii} = squeeze(imdata(rowsndx,colsndx,:));
				sumocolsX.images{ii} = squeeze(sum(imdata(rowsndx,colsndx,:),2,'omitnan'))';	
%				sumorowsY.images{ii} = squeeze(fixNaNimage(imdata(rowsndx,colsndx,:)));
%				sumocolsX.images{ii} = squeeze(sum(fixNaNimage(imdata(rowsndx,colsndx,:),2)))';			
%				sumorowsY.images{ii} = squeeze(fixNaNimage(imdata(rowsndx,colsndx,:)));
%				sumocolsX.images{ii} = squeeze(sum(fixNaNimage(imdata(rowsndx,colsndx,:),2)))';
			
			elseif NX==1;  % this won't execute if NY==1, so will mess up if ROI is a point!
				sumorowsY.images{ii} = squeeze(sum(imdata(rowsndx,colsndx,:),1,'omitnan'))';
				sumocolsX.images{ii} = squeeze(imdata(rowsndx,colsndx,:));	
%				sumorowsY.images{ii} = squeeze(sum(fixNaNimage(imdata(rowsndx,colsndx,:),1)))';
%				sumocolsX.images{ii} = squeeze(fixNaNimage(imdata(rowsndx,colsndx,:)));	
%				sumorowsY.images{ii} = squeeze(sum(fixNaNimage(imdata(rowsndx,colsndx,:),1)))';
%				sumocolsX.images{ii} = squeeze(fixNaNimage(imdata(rowsndx,colsndx,:)));	
	
			end
		
		else  %NY=NX=1 or a point ROI so it hardly makes any sense to sum!
		
				sumorowsY.images{ii} = squeeze(imdata(rowsndx,colsndx,:))';
				sumocolsX.images{ii} = sumorowsY.images{ii};
				
		end


end
end
%%%%%%%%%%%% HELPER FUNCTIONS AFTER THIS LINE %%%%%%%%%%%%
% function [Z]=helper1(inputts)
%
% end



