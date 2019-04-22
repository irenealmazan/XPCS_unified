function [ImageNew,BadPixDocument,MarkImage]	= apply_badpixreplace(ImageMatrix,badpixtxt,FLAGS)
%[ImageNew,BadPixDocument,MarkImage]	= apply_badpixreplace(ImageMatrix,badpix_filename,FLAGS)
% Used bad pix file ([Xbad Ybad Xrep Yrep] epics indexing of image matrix) 
%  to replace bad pixels with others in the image matrix 
%		Assume convention in badpix file is using (epic/imageJ) indexing
%			 enumerated as [X,Y] (start at 0)
%
%				(Matlab indexing would take these numbers   [Y+1,X+1]  (start at 1))
%
%			Matlab uses matrix indexing, so [Y=ROW, X=COLUMN]
%			And pixel [0,0] imageJ is matrix entry [1,1] in matlab
% 
%	OUTPUT
%		ImageNew  takes ImageMatrix and uses a simple pixel replacement as indicated by the badpix file
%		MarkImage can be used to indicate where bad pixels via an RGB indexing 
%			MarkImage(:,:,1) 1(bad pix) 0(good pix)
%			MarkImage(:,:,2) 1(pix used for replacement) 0(other)
%			MarkImage(:,:,3) zeros; RGB indexing >> image(MarkImage)
%			>> image(MarkImage) make figure on black ground showing
%				Red  pixel - are listed as BAD in pixel file
%				Green pixel  - are used for replacement to a bad pixel
%				Yellow pixel - are listed as BAD and used as a replacement
%				(that is, something has been typed in wrong in pixel file)
%	INPUT
%		ImageMatrix (can be 3 dimensional sequence of images)
%		badpix_filename  (string) e.g., 'badpix.txt' from EPICS or the matrix itself
%			each row is 4 numbers, [Xbad Ybad Xrep Yrep] 
%			e.g., [5  10  5  11] use pixel one column over as replacement
%				If using the NaN replacement option, it will ignore last columns but assume as Npoints x 2 array [Xbad1 Ybad1;Xbad2 Ybad2;...]
%		FLAGS (optional vector) default is [0 0]
%				FLAGS(1) 0= noplot(default) (1 plot image(MarkImage)) (useful to show where the badpix are in file)
%				FLAGS(2) 0  use the replacement (the way epics/imageJ would use the badpix file)
%						FLAGS(2)=1  instead replace badpix with NaN values
%		
if nargin<3;FLAGS=[0 0];end
if numel(FLAGS)==1;FLAGS(2)=0;end   % for legacy
ImageJ = 1 % assume that input is in spec/image format

[NYrow,NXcol,NZ] = size(ImageMatrix) %**

if isstring(badpixtxt);
	BADPIXTXT = load(badpixtxt);
else
	BADPIXTXT = badpixtxt;
end	
	[NUMbad,DUM]=size(BADPIXTXT);

% Find the indices for the bad pixels (note use in matlab, need to +1 indexing from what is in the badpix file
% Turns out that to index, cannot simply use A(Brow,Ccol) with Brow and Ccol vectors as it makes all combinations, use sub2ind instead
% Note that the '3' is used for consistency with rgb indexing used later, but is basically ignored 
	BADPIX	= BADPIXTXT(:,[1 2])+ImageJ;
	N1 		= sub2ind([NYrow NXcol 3],BADPIX(:,2),BADPIX(:,1));
	
	% Turns out that to index, cannot do A(Brow,Ccol) with Brow and Ccol vectors
	%	it makes all combinations between the two
	%	so turn the indices, [Brows,Ccol] into the linear index used in (A(:))
	N1 = sub2ind([NYrow NXcol 3],BADPIX(:,2),BADPIX(:,1),ones(NUMbad,1));
N2 = sub2ind([NYrow NXcol 3],REPPIX(:,2),REPPIX(:,1),ones(NUMbad,1).*2);
	
	% create the RGB indexed color map used for plotting the good/bad/ overlapped pixels
	% We do this even if we do not plot as it will also indicate when a pixel is listed as BAD as well as a REPLACEMENT
	%			(which happens when the badpix files pixels start accumulating without attention)
	MarkImage 			= zeros(NYrow,NXcol,3); % the '3' is used for RGB indexing
	MarkImage(N1)		= ones(NUMbad,1);		% Mark 'bad' pixels red (:,:,1) slice
	MarkImage(N2)		= ones(NUMbad,1);   	% Mark 'replace' pixels' green (N1.*2) puts the same indices into the (:,:,2) slice of the RGB
%	MarkImage(:,:,3) = zeros(NYrow,NXcol);		% Unnecessary, already zero 

if FLAGS(2)~= 1

	REPPIX = BADPIXTXT(:,[3 4])+ImageJ;   % the replacement pixels
	N2 = sub2ind([NYrow NXcol],REPPIX(:,2),REPPIX(:,1));
	N2 = sub2ind([NYrow NXcol 3],REPPIX(:,2),REPPIX(:,1),ones(NUMbad,1).*2);
	
	% make the replacement
	ImageNew = ImageMatrix;	
	for ii=1:NZ;
		Imageii = ImageMatrix(:,:,ii);
		Imageii(N1) = Imageii(N2);   
		ImageNew(:,:,ii) = Imageii;
	end
	
	BadPixDocument = char(['Replace badX badY with goodX goodY'],int2str(BADPIXTXT));

	% however, if doing a replacement, then we should check and tell the user the following if it happens!
	Nwrong = find([MarkImage(:,:,1)+MarkImage(:,:,2)]==2);   % both red and green!  (will be yellow on plot)
	if ~isempty(Nwrong);
		disp(['There are ',int2str(length(Nwrong)),' pixels that are listed both as BAD'])
		disp('and as good REPLACEMENT for another. That cannot be right')
		disp('Check and fix logic in the the bad pixel file')
	end


	
else   %simply replace with NaN for all NZ
		ImageNew 		= reshape(ImageMatrix,[],NZ);
		ImageNew(N1,:)	= NaN;
		ImageNew		= reshape(ImageNew,NYrow,NXcol,NZ);
	BadPixDocument = char(['NaN inserted in badX badY pixels'],int2str(BADPIXTXT(:,1:2)));
end	

% plot the helper image to show bad pix and replacement pixes
% Since it is easy to see NaN pixels with a simple imagesc(ImageMatrix(:,:,1), we don't process the plot option for that case
if FLAGS(1)==1 %  if all[FLAGS(1)==1 FLAGS(2)==0]
	figure;  N=gcf;
	H = image([1:NXcol]-1,[1:NYrow]-1,MarkImage);   
 % can plot each one separately, but together can be  RGB also
	xlabel('X Pixels Columns (in imageJ indexing)');
	ylabel('Y Pixels Rows (in imageJ indexing)');
	title(strvcat(	['Badpix file used [',pfilename(badpixtxt),']'],...
				'[RED  listed as BAD pixel]',...
				'[GREEN listed to replace some bad pixel]',...
				'[YELLOW listed as BAD but marked GOOD to replace another!]') )
	%setprops(H,NYrow,NXcol);
end

end	


