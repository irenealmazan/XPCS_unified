function FlatField = make_flatfield(specfilename,SCNNUM,loopscanFLAG,EXTRASTUFF);
%   older (2017_04) called as make_flatfield(str1,numsec,numpts,specfilename,SCNNUM);
%% str realted to equalization, numsec exposure, numpts (vector? of points) also filenameroot
%		now, put the extra stuff at end assuming that we have fixed sopme issues
%		not implemented at present in this version in order to get 
% used in 2017-07 beamtime using medipix 
% see book XXX p XXX
% since this was done in loopscan, could go back to the normal reading of files
% It can be collected 3 ways (depending on how crude we have the connections
% A) Taken with a loopscan, each exposure with each point in scan
%	loopscanFLAG=1; specFLAG=1;
% B) Taken with a loopscan, but on its own timing not exactly coincident with spec
%	loopscanFLAG=1; specFLAG=0;     
% C) Taken without a loopscan, without monitoring the monitor etc
%	loopscanFLAG=0;specFLAG=0
%%%  2017_0809_2 #64'
MOTORS = char('s6vgap','s6hgap','TY','TX','TY');
if nargin<2;disp('read the file to figure out what is needed');end
if nargin<3;loopscanFLAG = 1;end;

 [SPECpath,IMAGEpath,COMMONpath,HOMEpath] = pathdisplay;

% first 2 related to B and C not implemented right now
if ~loopscanFLAG;  % would be fixing up in here, likely need another flag passed to cover both situations
	specFLAG=0;
	SPECpath=[];
	MONaverage = 1;
	POSITIONS = NaN .* ones(length(MOTORS(:,1)),1);
	% Hardwired form of the file used book 197 page 75 for flatfield
	filenameroot = ['Th_' str1, '_' int2str(numsec) 'sec'];

	% had to specify manually if using the MPX3 stuff
	FlatField.filenameroot	= filenameroot;
	FlatField.seconds		= numsec;
	FlatField.Npts			= numpts;
	FlatField.monitor		= MONaverage;
	
	addnames2matrix
	[NameMatrix,sdata] = make_imnames_MPX3EPICS(filenameroot,numpts,NDIGITS,EPICStifpath);

elseif 0 % placeholder for the flag and timed coincindence 9will have to look back into the 2017_04 version
	specFLAG=1;
	fullspecfilename = [SPECpath,filesep,specfilename];
	[DATA,Npts,SCNTYPE,SCNDATE,ncols,collabels,fileheader,scnheader,comments] ...
				= readspecscan(fullspecfilename,SCNNUM);
	[OM,scninfo,OMspecial] = readspec_sevc_v6_helper(scnheader,fileheader,comments);
	INFOndx = chan2col(scninfo.allangles_label,MOTORS);
	POSITIONS = scninfo.allangles_i(INFOndx);
	MONaverage = mean(DATA(:,chan2col(collabels,'hexmon')));
	numsec = mean(DATA(:,chan2col(collabels,'Seconds')));
	T = table(MOTORS,POSITIONS);  disp(T) % prints out on screen
	
		% had to specify manually if using the MPX3 stuff
	FlatField.filenameroot	= filenameroot;
	FlatField.seconds		= numsec;
	FlatField.Npts			= numpts;
	FlatField.monitor		= MONaverage;
	

	[NameMatrix,sdata]=make_imnames_2017_04(filename,SCNNUM);
	
elseif loopscanFLAG;  % used currently since we had it with flag and timed coincinden 
	specFLAG=1;
	% this is a bit duplicative, sdata has a lot out, but not the 
	[NameMatrix,sdata]=make_imnames_2017_04(specfilename,SCNNUM);
		[OM,scninfo,OMspecial] = readspec_sevc_v6_helper(sdata.scnheader,sdata.fileheader,sdata.comments);
		INFOndx = chan2col(scninfo.allangles_label,MOTORS);
		POSITIONS = scninfo.allangles_i(INFOndx);
		T = table(MOTORS,POSITIONS);  disp(T) % prints out on screen
	MONITOR = sdata.DATA(:,chan2col(sdata.collabels,'hexmon'));  %%%%
	SECs    = sdata.DATA(:,chan2col(sdata.collabels,'Seconds'));  %%%%
	% in some sense, the next is legacy from 2017-04 way of doing it when we did not have loop scan coincindent with each image
	MONaverage = mean(MONITOR);
	numsec = mean(SECs);
	
	FlatField.filenameroot	= [NameMatrix.inputs.filename,' #',int2str(NameMatrix.inputs.scnnum)];;
	FlatField.seconds		= SECs;
	FlatField.Npts			= sdata.Npts;
	FlatField.monitor		= MONITOR;
end

%document some stuff
if specFLAG
	MONstr = ['<hexmon>ave=' num2str(MONaverage) ' cts/' int2str(numsec) 'sec'];
	SLITstr = ['s6vgap = ' num2str(prec(POSITIONS(1))) 'mm : s6hgap ' num2str(prec(POSITIONS(2))) 'mm'];
	SPECinfo = char([(specfilename) ' #' int2str(SCNNUM) ' : ' sdata.SCNTYPE],MONstr,SLITstr);
else
	SPECinfo = ['no spec file read'];
end
	
% read the files   % have to put put in the path	

[ImageStack,imtime] = load_MPX3(addnames2matrix([IMAGEpath,filesep],NameMatrix.fullfilenames));

% total (normalize to per image, if spec file input, also normalize to the tiff
sumImageStack = sum(ImageStack,3)./(FlatField.Npts);

% real would be 4096 or less - there are points that are a lot higher (1e5, for example)
% this takes the 'mean' of only those points that are 4096 or less
% so doesn't use the bad pixels that are 'high'
ImageMean = sum( [sumImageStack(:)<4096] .* sumImageStack(:) ) ./ sum(sumImageStack(:)<4096);

FlatField.MonMean		= MONaverage;
FlatField.ImageMean		= ImageMean;

if 0  %--- this was for the ones collected via MPX type programs
FlatField.Docu			= char(...
	[str1 ' : Equalization type'],...
	[int2str(numsec) ' sec : each exposure'],...
	['[' int2str(min(NameMatrix.points)) ':' int2str(max(NameMatrix.points)) '] : Points included'],...  %% this points field is not in the other
	['Root EPICS tif name is : ',FlatField.filenameroot],...
	['Normalization is to total number of tiffs : ' int2str(length(numpts))],...
	['ImageMean is ',num2str(ImageMean) ' cts/pix using only pix < 4096 cts'],...
	[SPECinfo]);
end
if 1 %--- this was for the ones collected 'normally'
FlatField.Docu		= char(...
	['? : Equalization type'],...
	[int2str(numsec) ' sec : each exposure'],...
	['File and scan in spec associated : ',FlatField.filenameroot],...
	['Normalization is to total number of tiffs : ' int2str(FlatField.Npts)],...
	['ImageMean is ',num2str(FlatField.ImageMean) ' cts/pix using only pix < 4096 cts'],...
	[SPECinfo]);
end


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
		
if 0  % old one
		TITLE  = char(...
					[pfilename([filenameroot '_  ']) int2str(min(Filenames.points)) ' to ' int2str(max(Filenames.points)) ': Flatfield'],...
					[Intensitystr],...
					pfilename(SPECinfo));
end
		TITLE  = char(...
					['Flatfield: ',Intensitystr],...
					[pfilename(SPECinfo)]);
					
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
	

	



	
