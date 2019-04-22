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
%					(note, there is option here to replace with NaN and not with neighboring pixel)
%					If the NaN is used, columns 3 and 4 (if they exist) will be ignored
%				If using the NaN replacement option, it will ignore last columns but assume as Npoints x 2 array [Xbad1 Ybad1;Xbad2 Ybad2;...]
%		FLAGS (optional vector) default is [0 0]
%				FLAGS(1) 0 noplot(default) (1 plot image(MarkImage)) (useful to show where the badpix are in file)
%				FLAGS(2) 0 replace badpix with NaN values (default);  =1 use the replacement (the way epics/imageJ would use the badpix file)
%							use the replacement (the way epics/imageJ would use the badpix file)
	
% eventually, a badpix file may not have the 'replacement' columns so assume 
% that we always have the badpix columns (column 1 and 2) and the replace columns (3 and 4) might not be there
%		(unless we are doing the replace pixels option as it done at the beamline)

% default, do not show the plot with the pixels marked, and replace Bad Pix with NaN (don't do the 'replace' with neighboring pixels
if nargin<3;FLAGS = [0 0];end
ImageJ = 1; 	% assume that the badpix is indexed as spec or imageJ (start at 0)

	[NYrow,NXcol,NZ]	= size(ImageMatrix);
	ImageNew			= ImageMatrix;

	BADPIXTXT			= load(badpixtxt);
	[NUMbad,Ncolumns]	= size(BADPIXTXT);

	% matrix (matlab) indexing starts at 1, 
	% ROI and imageJ used for badpix file start at 0 (assume always true)
	BADPIX			= BADPIXTXT(:,[1 2]) + ImageJ;
	REPPIX			= BADPIXTXT(:,[3 4]) + ImageJ;

	% Turns out that to index, cannot do A(Brow,Ccol) with Brow and Ccol vectors
	%	it makes all combinations between the two
	%	so turn the indices, [Brows,Ccol] into the linear index used in (A(:))
	N1		= sub2ind([NYrow NXcol 3],BADPIX(:,2),BADPIX(:,1),ones(NUMbad,1));
	% N2 = sub2ind([NYrow NXcol 3],REPPIX(:,2),REPPIX(:,1),ones(NUMbad,1).*2);

	% Create the MxNx[3] RGB indexed file used for visualization
	MarkImage		= zeros(NYrow,NXcol,3);
	MarkImage(N1)	= ones(NUMbad,1); %mark bad pixels as 'red'


ImageNew = ImageMatrix;

% Use replacement pixels  in the badpix file
if all([FLAGS(2)==1 Ncolumns==2]);
	
	% locate where to find the replacment pixels
	N2 = sub2ind([NYrow NXcol 3],REPPIX(:,2),REPPIX(:,1),ones(NUMbad,1).*2);
	
		% Additional information in visualization image
		MarkImage(N2) = ones(NUMbad,1);  % mark replacement pixels as 'green'

		% These will show up as 'yellow' in the visualization (r+g=yellow)
		Nwrong = find([MarkImage(:,:,1) + MarkImage(:,:,2)]==2);
		if ~isempty(Nwrong);
			disp(['There are ',int2str(length(Nwrong)),' pixels that are listed both as BAD'])
			disp('and as good REPLACEMENT for another. That cannot be right')
			disp('Check and fix logic in the the bad pixel file')
		end

	% perform the Replacements (using the replacement pixels indicated in badpixfile
	for ii = 1:NZ;
		Imageii 		= ImageMatrix(:,:,ii);
		Imageii(N1) 	= Imageii([N2 - (NYrow.*NXcol)]);
		ImageNew(:,:,ii) = Imageii;
	end
	
	VISstr = char(['Badpix file used [', pfilename(badpixtxt),']'],...
				'[RED  listed BAD pixel]',...
				'[GREEN listed to replace some bad pixel]',...
				'[YELLOW Hey - BAD and marked as good to replace another]');
				
	BadPixDocument	= char(['Replace badX badY with goodX goodY  points below '],int2str(BADPIXTXT));
	DOCU			= char(['Replace badX badY with goodX goodY using info from '  (badpixtxt)]);
	disp(DOCU);

else   % put NaN in the bad pixel positions


	% perform the Replacements (using the replacement pixels indicated in badpixfile
	for ii = 1:NZ;
		Imageii 		 = ImageMatrix(:,:,ii);
		Imageii(N1) 	 = NaN .* ones(size(N1));
		ImageNew(:,:,ii) = Imageii;
	end

	VISstr = char(['Badpix file used [', pfilename(badpixtxt),']'],...
				'[RED  listed BAD pixel]');
				
	BadPixDocument 	= char(['Replace badX badY with NaN using points below '],int2str(BADPIXTXT));
	DOCU			= char(['Replace badX badY with NaN using info from ' (badpixtxt)]);
	disp(DOCU);

end

% If requested by FLAGS(2); show the visualization of the bad pixel positions
% Visualization will be 'black' with red for bad pixels
%		(and if replacements, then green for good, and yellow for ones that are listed as bad and as a replacement for another pixel

if ~(FLAGS(2)==1);
	figure;  N1=gcf;
		H1 = image([1:NXcol]-ImageJ,[1:NYrow]-ImageJ,MarkImage);   
		% can plot each one separately, but together can be  RGB also
		xlabel('X Pixels Columns');
		ylabel('Y Pixels Rows');
		title(VISstr);
%		plot_adjust(gca,0);
	%setprops(H,NYrow,NXcol);
	%subplot(3,1,2);
	figure; N2=gcf;
		H2 = imagesc([1:NXcol]-ImageJ,[1:NYrow]-ImageJ,sum(ImageMatrix,3));   
		A2 = gca;
%		plot_adjust(A2,0);
		% can plot each one separately, but together can be  RGB also
		xlabel('X Pixels Columns');
		ylabel('Y Pixels Rows');
		VISstr2 = char(['The Image data without badpixel replacements'],...
				'(summed over 3rd dimension if it is sequence of images');
		title(VISstr2);
	%subplot(3,1,3);
	figure; N3=gcf;
		H3 = imagesc([1:NXcol]-ImageJ,[1:NYrow]-ImageJ,sum(ImageNew,3)); 
%		MM = mean(mean(mean(ImageNew)));  
		% can plot each one separately, but together can be  RGB also
%		plot_adjust(gca,0);
		xlabel('X Pixels Columns');
		ylabel('Y Pixels Rows');
		VISstr3 = char(['To help check for bad, summing the Image data entered'],...
				pfilename(DOCU));
		title(VISstr3)
		CC = get(gca,'clim');
		set(A2,'clim',CC);
		

end

end	  % of function


