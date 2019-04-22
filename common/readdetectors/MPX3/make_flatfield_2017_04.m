function FlatField = make_flatfield(str1,numsec,numpts,specfilename,SCNNUM);
% used in 2017-07 beamtime using medipix 
% see book 197 page 92

% hardwired
NDIGITS = 6;  % 000000 000001 etc
EPICStifpath = '/home/cthompson/DATA/2017/2017_04_mocvd_nitrideXPCS/areadata/medipix_3_EPICS/EPICS_next';

MOTORS = char('s6vgap','s6hgap','TY','TX','TY');
if nargin< 5;
	specFLAG=0;
	SPECpath=[];
	MONaverage = 1;
	POSITIONS = NaN .* ones(length(MOTORS(:,1)),1)
else 
	specFLAG=1;
	SPECpath = pathdisplay;
	fullspecfilename = [SPECpath,filesep,specfilename];
	[DATA,Npts,SCNTYPE,SCNDATE,ncols,collabels,fileheader,scnheader,comments] ...
				= readspecscan(fullspecfilename,SCNNUM);
	[OM,scninfo,OMspecial] = readspec_sevc_v6_helper(scnheader,fileheader,comments);
	INFOndx = chan2col(scninfo.allangles_label,MOTORS);
	POSITIONS = scninfo.allangles_i(INFOndx);
	MONaverage = mean(DATA(:,chan2col(collabels,'hexmon')));
	SECspec = mean(DATA(:,chan2col(collabels,'Seconds')));
	T = table(MOTORS,POSITIONS);  disp(T) % prints out on screen
end

%document some stuff
if specFLAG
	MONstr = ['<hexmon>=' num2str(MONaverage) ' cts/' int2str(SECspec) 'sec'];
	SLITstr = ['s6vgap = ' num2str(prec(POSITIONS(1))) 'mm : s6hgap ' num2str(prec(POSITIONS(2))) 'mm'];
	SPECinfo = char([(specfilename) ' #' int2str(SCNNUM) ' : ' SCNTYPE],MONstr,SLITstr);
else
	SPECinfo = ['no spec file read'];
end
	
% Hardwired form of the file used book 197 page 75 for flatfield
filenameroot = ['Th_' str1, '_' int2str(numsec) 'sec'];
Filenames = make_imnames_MPX3EPICS(filenameroot,numpts,NDIGITS,EPICStifpath);
	
[ImageStack,imtime] = load_MPX3(Filenames.fullfilenames);

% total (normalize to per image, if spec file input, also normalize to the tiff
sumImageStack = sum(ImageStack,3)./length(numpts);

% real would be 4096 or less - there are points that are a lot higher (1e5, for example)
% this takes the 'mean' of only those points that are 4096 or less
% so doesn't use the bad pixels that are 'high'
ImageMean = sum( [sumImageStack(:)<4096] .* sumImageStack(:) ) ./ sum(sumImageStack(:)<4096);

FlatField.filenameroot	= filenameroot;
FlatField.seconds		= numsec;
FlatField.Npts			= numpts;
FlatField.MonMean		= MONaverage;
FlatField.ImageMean		= ImageMean;
FlatField.Docu			= char(...
	[str1 ' : Equalization type'],...
	[int2str(numsec) ' sec : each exposure'],...
	['[' int2str(min(Filenames.points)) ':' int2str(max(Filenames.points)) '] : Points included'],...
	['Root EPICS tif name is : ',filenameroot],...
	['Normalization is to total number of tiffs : ' int2str(length(numpts))],...
	['ImageMean is ',num2str(ImageMean) ' cts/pix using only pix < 4096 cts'],...
	[SPECinfo]);

FlatField.imtime 		= imtime;

FlatField.flatfield = sumImageStack;
FlatField.flatfieldscaled = FlatField.flatfield ./ FlatField.ImageMean;

XCOLpts = [0:515];YROWpts=[0:515];
AXISdet = [250 510 250 510];
AXISdet = [0 515 0 515];
		
		
		figure;clf
		set(gcf,'Name','Summed Flatfield');
		
		
%		Intensitystr = 'Log(Summed(I)/points)';
%		HS = pcolor(XCOLpts,YROWpts,log(sumImageStack(:,:)) );

		Intensitystr = ['I=Isum/(Mean * Npts) : <Imean(<4096)>=' num2str(ImageMean) ' cts/pix'];
		HS = surf(XCOLpts,YROWpts,(FlatField.flatfieldscaled(:,:)) );
		colorbar
			if ~isempty(AXISdet); axis(AXISdet);end
			%if ~isempty(CLIMlin); 
			%	set(gca,'clim',CLIMlin);
				%DOCUclim = mDOCUclim(CLIMlin);
			%end;
		
		shading flat
		axis square
		view(0,90);
		set(gca,'clim',[0 2]);
		
		DOCUclim = mDOCUclim(gca);
		
		TITLE  = char(...
					[pfilename([filenameroot '_  ']) int2str(min(Filenames.points)) ' to ' int2str(max(Filenames.points)) ': Flatfield'],...
					[Intensitystr],...
					pfilename(SPECinfo));
		title(TITLE);
		xlabel('columns (X)'),ylabel('rows (Y)');
		prettyplot(gca);
		set(gcf,'paperposition',[2 2 5 5])

end
%% helper functions %%%%%%%%%%%%%%%%%%%%%
function Xout = prec(Xin,DIG)
% change 0.010123333 to 0.0101 (rounding conventional)
	if nargin<2;DIG = 4;end
	
	Xout = round(Xin .* 10^DIG)./10^DIG;
	
end
	

	



	
