
[SPECpath,AREApath,COMMONpath,HOMEpath] = pathdisplay;
eval(RUNPARAMS);

%specfilename = '2016_1015_2';
%HDFfilename = specfilename;
%DOCU0  = 'n161012a (AH1385.10.2) m-plane GaN';

clear Iwin Iwinnorm;
HLine=[];HSurf=[];

if HDF5==1   % not really used, but legacy for when using pixirad with windows
	STR.scanflag=1;
	STR.imname=[];
	STR.SPECpath=DATApath;
	STR.AREApath=AREApath;
	STR.p_image=[];
	STR.ending='.h5';
else  % medipix, pilatus, all saved in tif thus far (usually)
	STR.scanflag=1;
	STR.imname=[];      % image names use standard unique naming schemes
%	STR.imname='pil_';  % caxis non-unique naming scheme for images (sigh) 
	STR.SPECpath=DATApath;
	STR.AREApath=AREApath;   % need to change pathdisplay if need pilatus
	STR.p_image=['p_image'];   %% [] if it uses the preferred points instead
	STR.ending='.tif';
end

ImageJ = 1; % true - using imageJ indexing [start at 0] when denoting ROIs and when plotting pixels
%	This applies to AXISdet and ROI values so must be consistent.

% FORCE only one scan at a time
SCNs = eval(SCNstr);SCNs = SCNs(1); 

if strncmp('caxis',SCNFUNC,min([length('caxis') length(SCNFUNC)]));  % which diffractometer
	readspechelp = str2func('readspec_caxis_v4_helper');
	INFOinterest = ['Sam_Z  stemp  Chi  Eta'];     %!!! put 2 spaces between each!!!!
else
	readspechelp = str2func('readspec_sevc_v6_helper');
	INFOinterest = ['s6vgap  s6hgap  Eta  Phi  Chi  Nu  Delta  TY  TX  TZ'];     %!!! put 2 spaces between each!!!!
end

%----------  Make the matrix of the names for the images, note sdata holds information from spec file
[NameMatrix,sdata]	= make_imnames_2017_07(HDFfilename,SCNs,STR);
FullNameHDF 		= addnames2matrix([STR.AREApath,filesep],NameMatrix.fullfilenames);

%----------  At this point can pull out quite a bit from the spec file (in sdata)
	[OM,scninfo,OMspecial]	= readspechelp(sdata.scnheader,sdata.fileheader,sdata.comments);
	INFOndx 				= chan2col(scninfo.allangles_label,INFOinterest);
	INFOpositions 			= scninfo.allangles_i(INFOndx);
		DUM=[];
		% put into a line
		for ii=1:length(INFOpositions);
			DUM = [DUM,pdigit(INFOpositions(ii)) ' : '];
		end
	INFOstr = char(INFOinterest,DUM);   % currently added onto 0 and 5
	% INFOstr = [];    % turn this off quickly if desired 
	
	secperpoint		= sdata.DATA(:,chan2col(sdata.collabels,timename));
	timestampSpec	= sdata.DATA(:,chan2col(sdata.collabels,etimename));
	timestampEpoch  = sdata.DATA(:,chan2col(sdata.collabels,epochtimename));
	filterpoint  	= sdata.DATA(:,chan2col(sdata.collabels,filtername));
	Monitor			= sdata.DATA(:,chan2col(sdata.collabels,monname));
	MONave 			= mean(Monitor(:,1));  %% occasionally counters in more than one column
	mon100mA 		= MONave./mean(secperpoint);  %do not use RUNPARAM- since slit size changes a lot

	%---- get the normalization corrections ready
	if NORMFLAG == 1;	
	Norm = MONave .* filter_correction(filterpoint,filter1,filtercorr)./(secperpoint.*Monitor(:,1));
	NORMdoc = ['Norm = monave*filtercorrection/(sec per point * Monitoer per point): WITH filter correction']
	else 
	Norm = MONave ./(secperpoint.*Monitor(:,1));
	NORMdoc = ['Norm = monave/(sec per point * Monitoer per point): no filter correction']
	end



%---------   Read the images
if HDF5~=0;  % assume HDF5 was pixirad in windowing form
	[Iwin,timestampX,IColors] = loadhdf5_CT(FullNameHDF,[LOWHIGH]);
	timestampX_flag = 1;  % time stamp flag has good time resolution and is precise
	% also cannot be changed by resaving file
	THRESHusedoc = [' {',THRESHuse, '} '];
else
	[Iwin,timestampX] = load_MPX3(FullNameHDF);
	% The following are not used for tifs, only for pixirad in window mode
	% Note the timestamp for the medipix is only to 1 second, not to milliseconds
	% So need to use the spec information for time.
	% timestamp for medipix also weird in 2017-08
	Icolors = [];
	THRESHuse = [];   %% hardwiring to not use the pixirad high low, etc
	THRESHusedoc = [];
	timestampX_flag = 0;   % timestamp flag from tif is not great, keep using spec
end


if ~LOGFLAG
	DOCUInt = ['Inorm',THRESHusedoc];
else
	DOCUInt = ['log10(Inorm',THRESHusedoc,')'];
end

if FLATFIELDFLAG
	DOCUInt = [DOCUInt,'[FF]'];
end

TITLEstr1 = char(...
		[pfilename(HDFfilename),' #', SCNstr,' : ', sdata.SCNDATE],...
		[sdata.SCNTYPE],...
		[pfilename(DOCU0),' ',DOCUscan]);

TITLEstrshort = char(...
		[pfilename(HDFfilename),' #',SCNstr,' : ',DOCUInt],...
		[sdata.SCNTYPE]);
		

%------------   If keeping track of valves to plot
if ~strncmp('caxis',SCNFUNC,min([length('caxis') length(SCNFUNC)]));  % which diffractometer

	valve_p	= sdata.DATA(:,chan2col(sdata.collabels,'valve_p'));
	valve_v	= sdata.DATA(:,chan2col(sdata.collabels,'valve_v'));
	lastframes = find(diff(valve_p)~=0);
	disp(['Last frames (Spec point number convention) before valve switches: ' num2str(lastframes')]);
end  % NEED TO FIX LATER SO ALSO WORKS WITH CAXIS
	
		if timestampX_flag==0;
			timestampX = timestampSpec + timestampEpoch(1);
		end
	

%	if isempty(Xstepcol); Xstepcol=1;end   % use scan variable
%	 
%		if Xstepcol<1;
%			SCNXLABEL = 'spec point # ';
%		else
%			SCNXLABEL	= col2chan(sdata.collabels,Xstepcol);  
%		end
%	
%		if strncmp(SCNXLABEL,'Time',4); 
%			SCNXLABEL	= [SCNXLABEL,' [\Delta sec] from SpecFile '];  % in pixirad, 
%			Xsteps	= timestampX-timestampX(1);
%		elseif Xstepcol<1;
%			Xsteps	= [0:length(timestampX)-ImageJ];  
%		else
%			Xsteps	= sdata.DATA(:,Xstepcol);
%		end

	if FORCE.FLAG; %Xstepcol=1;end   % use scan variable	 
		if strfind(FORCE.X,'specpoint')
			SCNXLABEL = 'spec point # ';
			Xstepcol=-1;
		else
			SCNXLABEL	= FORCE.X;
			Xstepcol = chan2col(sdata.collabels,FORCE.X); 
				if isempty(Xstepcol);disp('ERROR: asking for Xaxis not present in file');end
		end	
	else
		Xstepcol = 1;
		SCNXLABEL = col2chan(sdata.collabels,1); 
	end
	
	if strncmp(SCNXLABEL,'Time',4); 
		SCNXLABEL	= [SCNXLABEL,' [\Delta sec] from SpecFile '];  % in pixirad, 
		Xsteps	= timestampX-timestampX(1);
	elseif Xstepcol<1;
		Xsteps	= [0:length(timestampX)-ImageJ];  
	else
		Xsteps	= sdata.DATA(:,Xstepcol);
	end

	

%	Here we are being consistent with enumeration in SPEC points and ImageJ
%		ImageJ convention is used in the SPEC ROI's and EPICS ROIs images
%		SPEC convention is used in its points when we 'write' down a data point
% ImageJ = 1 % using ImageJ convention, ImageJ=0 use matlab convention
%		matlab convention numbering starts from 1, not 0 
%			from the screen	
	SPECpts = [1:(length(Xsteps))]-ImageJ;
	YROWpts = [1:(length(Iwin(:,1,1)))]-ImageJ;  
	XCOLpts = [1:(length(Iwin(1,:,1)))]-ImageJ;
	
% make some of the images to be ready for other work, some is legacy from pixirad windowing

		if strncmp(THRESHuse,'low',3)
			II = IColors.low;
		elseif strncmp(THRESHuse,'hig',3)
			II = IColors.high;
		else
			II = Iwin;
			%THRESHuse = 'win';  % default if pixirad in windowed versoin
			THRESHuse = [];
		end
				
		for ii=1:length(Iwin(1,1,:));
			IInorm(:,:,ii)		= II(:,:,ii) .* Norm(ii);
			Iwinnorm(:,:,ii)	= Iwin(:,:,ii) .* Norm(ii);
		end
		% Calculate mean background image
		if isempty(BKG)
			IIimbkg = 0;
			Iwinimbkg = 0;
		else
			IIimbkg		= mean(IInorm(:,:,[BKG(1):BKG(2)]+ImageJ),3);
			Iwinimbkg	= mean(Iwinnorm(:,:,[BKG(1):BKG(2)]+ImageJ),3);
		end
		
		% note if FLATFIELDflag=0; then user_input should make FLATFIELD=1;
		for ii=1:length(Iwin(1,1,:));		
			IInormb(:,:,ii)		= (IInorm(:,:,ii)	- IIimbkg)./FLATFIELD;
			Iwinnormb(:,:,ii)	= (Iwinnorm(:,:,ii)	- Iwinimbkg)./FLATFIELD;
		end



%%%%  CT took out the minus flags choices(used for the pixirad and windowing) %%%
%%%% and put in them in a separate script, I will comment it out as we are not using it
%%%% for the medipix anyway but it helps clean up this frankenstein file

%user_input_guts_MPX3_minusflags

%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(WhatToDo==0); 

	if ~LOGFLAG
		figure;clf
		set(gcf,'Name',['Lin Image',THRESHusedoc,' Spec # ' int2str(SINGLE)]);
		HS = pcolor(XCOLpts,YROWpts,(IInormb(:,:,SINGLE+ImageJ)));
		shading interp
		HSurf.IInormb = HS;

			if ~isempty(AXISdet); axis(AXISdet);end
			if ~isempty(CLIMlin); 
				set(gca,'clim',CLIMlin);
				DOCUclim = mDOCUclim(CLIMlin);
			end;
	else
		figure;clf
		set(gcf,'Name',['Log10 Image',THRESHusedoc,' Spec #' int2str(SINGLE)]);
		HS = pcolor(XCOLpts,YROWpts,log10(fixlogimage((IInormb(:,:,SINGLE+ImageJ)))));
		shading flat
		HSurf.IInormb = HS;

			if ~isempty(AXISdet), axis(AXISdet);end
			if ~isempty(CLIMlog); 
				set(gca,'clim',CLIMlog);
				DOCUclim = mDOCUclim(CLIMlog);
			end;
	end
	
	DOCUlast = ['(Spec pt #',int2str(SINGLE),') ', DOCUInt, DOCUclim];
	title(char(TITLEstr1,DOCUlast,pfilename(INFOstr)));
	xlabel(XCOL);ylabel(YROW);
	prettyplot(gca,[2 2 4 3.8]);
end	


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(WhatToDo==1) %& ~isempty(ROIS) %% open up make ROIS
	myROIS = makerois(4,gcf);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(WhatToDo==2)

	[HROI,COLORORDER] = showrois(ROIS,gcf);
		
		figure;clf;
		set(gcf,'Name','Enumeration of ROIs colors');
		if ~isempty(AXISdet); axis(AXISdet);end
		prettyplot(gca);    % put here or it clobbers the showROIS colororder 		
		showrois(ROIS,gcf,2,COLORORDER);
		title('Enumeration of the ROIs and their colors');
		if ~isempty(AXISdet); 
			axis(AXISdet);
		else
			axis([min(XCOLpts) max(XCOLpts) min(YROWpts) max(YROWpts)]);
		end

			LEGEND = [char([ones(length(HROI),1)*'ROI #']) int2str([1:length(HROI)]')];
			legend(LEGEND);
		
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(WhatToDo==3)  %% make movie  (frame must be 560x413?)
	figure;clf;figure(gcf);
	set(gcf,'Name',['MOVIE ' THRESHusedoc]);
	% set up initial size and image and climits  (Default SINGLE is 0 (spec point + ImageJ = matlab point)
	%	
	% Note that video get upset if size changes of image some of this is to 
	% try to make sure this doesn't happen but it still happens and I 
	% don't know how to fix at times except to keep trying again (magic)
	clear HS;
	HS = pcolor(XCOLpts,YROWpts,log10(fixlogimage(IInormb(:,:,SINGLE+ImageJ))));
	DOCUInt3 = [' log10(Inorm) '];
	shading flat
		HT = title(...
			char([TITLEstrshort],...
			['Spec Pt # [ 00000 ] : step =  [0000000]',mDOCUclim(CLIMlog)]));	

		if ~ASPECTflag
		axis equal;
		else
			HA = gca;
			%daspect(gca,ASPECT);   %this might not be in earlier version
			HA = gca;,
			set(HA,'DataAspectRatioMode','manual','DataAspectRatio',[ASPECT]);
			end
			
		
		xlabel(XCOL);ylabel(YROW);
		
		if MOVIEROIFLAG
			MAA = ROIS(1,:);
		else
			if ~isempty(AXISdet);
				MAA = AXISdet;
			else
				MAA = [XCOLpts(1) XCOLpts(end) YROWpts(1) YROWpts(end)];
			end
		end	

		set(gca,'xlim',MAA([1 2]),'ylim',MAA([3 4]));
		
		if ~isempty(CLIMlog); 
			set(gca,'clim',CLIMlog);
		end;
		prettyplot(gca);
			if exist('COLORORDER','var')%length(ROIS(:,1))<length(COLORORDER(:,1))
				showrois(ROIS,gcf,0.5,COLORORDER);
			else
				showrois(ROIS,gcf,0.5);
			end

	% Prepare the new file.
		if exist('M1.avi')==2; delete('M1.avi'); end  
		vidObj = VideoWriter('M1.avi');
		vidObj.FrameRate = 5;		% default 30 fps


		open(vidObj);


%	for ii=1:length(timestampX)  % CT
	for ii=1:length(Xsteps)
		PTstr = int2str(ii-1+10000);PTstr = PTstr(2:end);
		
		if ~LOGFLAG
			set(HS,'CData',((IInormb(:,:,ii))));
			if ~isempty(CLIMlin); set(gca,'clim',CLIMlin);end;
		else
			set(HS,'CData',(log10(fixlogimage(IInormb(:,:,ii)))));
			if ~isempty(CLIMlog); set(gca,'clim',CLIMlog);end;
		end		
		set(HT,'String',...
			char([TITLEstrshort],...
			['Spec Pt # [',PTstr,'] : step = [',num2str(Xsteps(ii)),']',mDOCUclim(gca)]));
		currFrame = getframe(gcf);
		writeVideo(vidObj,currFrame);
	end
	disp('*')
	%Mov_n160208a_intval=M;
	close(vidObj);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(WhatToDo==4) %% make summed images 
%	ImageJ = 1;  % this must be set and be consistent earlier! 
%	ImageJ indexing (start (0,0)) in ROI indexing
	 [SIZE1,SIZE2,SIZE3]=size(IInormb);  % in case whole scan not complete
	 [SoY,SoX] = slicesumrois(IInormb,ROIS,ImageJ);%(1,:));
%	Note - This does not 'normalize' the slices to the number of pixels summed
%			However, after using function, use SoY.image{i}./SoY.norm{i} 
%				for ROI {i} 

	for ii=1%:length(ROIS(:,1));

  if SoXFLAG
	figure;clf
	set(gcf,'Name',['SoX ROIs #' int2str(ii)]);

	 if ~LOGFLAG;
	 HS = pcolor(Xsteps,SoX.ndx{ii},((SoX.images{ii})));
	 else
	 HS = pcolor(Xsteps,SoX.ndx{ii},log10(fixlogimage(SoX.images{ii})));
	 end

	 HPC.SoX_ROIS(ii) = HS;
	 shading flat; 
	 xlabel(SCNXLABEL);ylabel(YROW);
	 title(char(TITLEstr1,[DOCUInt, ' summed over XCOL in ROI #' int2str(ii)]));
	 			if ~isempty(AXISdet); set(gca,'Ylim',AXISdet(3:4));end
%	 			if ~isempty(AXISdet); set(gca,'Ylim',ROIS(ii,[3:4]));end
	 prettyplot(gca)
	 	 makeyline(Xsteps(lastframes));
  end
	
  if SoYFLAG  % don't plot if SoY flag 0
	figure;clf
	set(gcf,'Name',['SoY ROIs #' int2str(ii)]');
	
	 if ~LOGFLAG;
	 HS = pcolor(Xsteps,SoY.ndx{ii},((SoY.images{ii})));
		
	 else
	 HS = pcolor(Xsteps,SoY.ndx{ii},log10(fixlogimage(SoY.images{ii})));
	 end

	 HPC.SoY_ROIS(ii) = HS;
	 shading flat; 
	 xlabel(SCNXLABEL);ylabel(XCOL);
	 title(char(TITLEstr1,[DOCUInt ' summed over YROW in ROI #' int2str(ii)]));
	 % Note Ylim in this context is on output plot not the since we put the time etc on X axis
%	 			if ~isempty(AXISdet); set(gca,'Ylim',AXISdet(1:2));end
	 			if ~isempty(AXISdet); set(gca,'Ylim',ROIS(ii,[1:2]));end
	 prettyplot(gca) 
	 	 makeyline(Xsteps(lastframes));
	 
	 end
  end

if  1  
	if SoXFLAG
	figure; clf
	set(gcf,'Name','SoX and Sum Points ');
	 HL = semilogy(SoX.ndx{1},sum(SoX.images{1}'));
		for jj = 2:length(SoX.images);
			HL(jj) = line(SoX.ndx{jj},sum(SoX.images{jj}'));
		end


	 xlabel(YROW),ylabel('Int (arb)'),
	 title(char(TITLEstr1,['summed over XCOLS in all ROI and ',int2str(SIZE3),' scan pts']));
	 prettyplot(gca);  % 
		for ii=1:length(ROIS(:,1));
		set(HL(ii),'Color',COLORORDER(ii,:),'LineWidth',2);
		end
		
		HLine.SoX_SPts = HL;clear HL;
	
	end
	if SoYFLAG	
	figure; clf
	set(gcf,'Name','SoY and Sum Points ');
	 HL = semilogy(SoY.ndx{1},sum(SoY.images{1}'));
		for jj = 2:length(SoY.images);
			HL(jj) = line(SoY.ndx{jj},sum(SoY.images{jj}'));
		end

	 xlabel(XCOL),ylabel('Int (arb)'),

	 title(char(TITLEstr1,['summed over YROWS in all ROI and ',int2str(SIZE3),' scan pts']));
	 prettyplot(gca);  % 
		for ii=1:length(ROIS(:,1));
		set(HL(ii),'Color',COLORORDER(ii,:),'LineWidth',2);
		end
	 
		HLine.SoY_SPts = HL;clear HL;
	end

	if ~isempty(POINTSUMS)   % sum over points sets in POINTSUMS
	
		if SoXFLAG
		figure; clf;
		set(gcf,'Name','SoX 1st ROI and Sum select points');
		for ii=1:length(POINTSUMS(:,1));
			Ni = [POINTSUMS(ii,1): POINTSUMS(ii,2)]+ImageJ;  %use as matlab
			HL(ii) = line(SoX.ndx{1},sum(SoX.images{1}(:,Ni)'));
		end
		
		
			set(gca,'Yscale','log')
			xlabel(YROW),ylabel('Int (arb)'),
			title(char(TITLEstr1,['summed over X in 1st ROI and between selected Spec Points']));
			prettyplot(gca);legend(addnames2matrix('sum between [', int2str(POINTSUMS),']'))
		HLine.SoX_SoSelectPt = HL;clear HL;
		end
		
		if SoYFLAG

		figure; clf;
		set(gcf,'Name','SoY 1st ROI and Sum selected points');
		for ii=1:length(POINTSUMS(:,1));
			Ni = [POINTSUMS(ii,1): POINTSUMS(ii,2)]+ImageJ;
			HL(ii) = line(SoY.ndx{1},sum(SoY.images{1}(:,Ni)'));
		end
			set(gca,'Yscale','log')
			xlabel(XCOL),ylabel('Int (arb)'),
			title(char(TITLEstr1,['summed over Y in 1st ROI and between selected Spec Points']));
			prettyplot(gca);legend(addnames2matrix('sum between [', int2str(POINTSUMS),']'));
		HLine.SoY_SoSelectPt = HL;clear HL;
		end
	end
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% note if the ROI has NaN in any of the pixels, this will mess up and give NaN
%%% need to fix sumrois to fix that
if any(WhatToDo==5)  %% make summed lines of ROI's on images
	figure;clf
	set(gcf,'Name','ROIs as counters');
	 
	 [sumROIS,sumROISnormed] = sumrois_wNaN(IInormb,ROIS,ImageJ);
%	 [sumROIS,sumROISnormed] = sumrois(IInormb,ROIS,ImageJ);
	 %	Isum	Array  : each column for each ROI. Length of vectors (rows) num of images
	 %HL = semilogy(Xsteps,sumROIS);
	 HL = plot(Xsteps,sumROISnormed); YLABEL = 'Inorm(summed over ROI)';
%	 HL = plot(Xsteps,log10(sumROISnormed)); YLABEL = 'log10[Inorm(summed over ROI)/ROIsize]';

% mocvd
	 makeyline(Xsteps(lastframes)); % need to turn these into colors so do before prettyplot
	 xlabel(SCNXLABEL);ylabel(YLABEL);
	 title(char(TITLEstr1,INFOstr,YLABEL));
	 prettyplot(gca)
	 
	 for ii=1:length(ROIS(:,1));
		set(HL(ii),'Color',COLORORDER(ii,:),'LineWidth',2);
	 end
	 
	 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(WhatToDo==6);   %to make it automatic for summs
	if isempty(POINTSUMS); POINTS=[1 length(Xsteps)]-ImageJ;
	else POINTS = POINTSUMS;
	end
	
for ii = 1:length(POINTS(:,1))
	Ni = [POINTS(ii,1): POINTS(ii,2)]+ImageJ;
	NP = length(Ni);

	
	figure;clf;
	set(gcf,'Name','Image summed over points/Num Points');	
	if ~LOGFLAG

		HS = surf(XCOLpts,YROWpts,sum(Iwinnorm(:,:,Ni),3)./(NP));

			if ~isempty(AXISdet); axis(AXISdet);end
			
			if ~isempty(CLIMlin); 
				set(gca,'clim',CLIMlin);
				DOCUclim = mDOCUclim(CLIMlin);
			else 
				CLIMtemp = get(gca,'clim');
				DOCUclim = mDOCUclim(CLIMtemp);
			end;
			DOCUInt6 = [DOCUInt, ' SUMMED over SPEC scan pts [',int2str(POINTS(ii,:)),']'];
	else
	
		HS = surf(XCOLpts,YROWpts,log10(fixlogimage(sum(IInormb(:,:,Ni),3)./(NP))));
			if ~isempty(AXISdet), axis(AXISdet);end
			if ~isempty(CLIMlog); 
				set(gca,'clim',CLIMlog);
				DOCUclim = mDOCUclim(CLIMlog);
			else 
				CLIMtemp = get(gca,'clim');
				DOCUclim = mDOCUclim(CLIMtemp);
			end;
			DOCUInt6 = [DOCUInt, ' SUMMED over SPEC scan pts) [',int2str(POINTS(ii,:)),'])'];
	end
	
	view(0,90);shading flat;
	
	if FLATFIELDFLAG
		DOCUInt = [DOCUInt,' [FF]'];
	end
	
	DOCUlast = char(DOCUInt6, INFOstr);
	title(char(TITLEstr1,DOCUlast));xlabel(XCOL);ylabel(YROW);
	prettyplot(gca,[2 2 4 3.8]);
			if length(ROIS(:,1))<length(COLORORDER(:,1))
				showrois(ROIS,gcf,0.5,COLORORDER);
			else
				showrois(ROIS,gcf,0.5);
			end
	axis square
	
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(WhatToDo==7); % DO NOT USE Calculate two-time correlations of pixels in ROIs over all frames
    
    disp('Calculating 2-time correlations, please wait.');
    
    Nt = length(Xsteps);
    Iroi = ones(0,Nt); % get all pixels in ROIs
    for ii = 1:size(ROIS,1);
        Iroi = [Iroi; reshape(IInormb(ROIS(ii,3):ROIS(ii,4),ROIS(ii,1):ROIS(ii,2),:),[],Nt)];
    end
    
    IM = mean(Iroi)';    
    IIM = NaN*ones(Nt,Nt);
    for ii = 1:Nt
        for jj = 1:ii
            IIM(ii,jj) = mean(Iroi(:,ii).*Iroi(:,jj));
            IIM(jj,ii) = IIM(ii,jj);
        end
    end
    IID = diag(IIM) - IM.^2;
    CC = (IIM - IM*IM')./sqrt(IID*IID'); % Normalized to make diagonal unity
    CA = IIM./(IM*IM') - 1; % Normalized to make field = contrast, diagonal ~kbar^-1
    
    
    figure; % plot of kbar vs frame
    set(gcf,'Name','Plot of kbar vs frame (scan units or point)');
    	%axes('Box','on');  Prettyplot does thi
    HL = line(Xsteps,IM);

		xlabel(SCNXLABEL);
		ylabel('Mean cts / pixel');
		   title(char(['kbar vs frame'],TITLEstr1));
		   	 	 makeyline(Xsteps(lastframes));
		prettyplot(gca);
    
    figure;
    set(gcf,'Name','2-time averaged all ROIS');
    HS = pcolor(CC);shading flat;
		if ~isempty(CLIMXP);
			set(gca,'clim',CLIMXP);
		else
			CLIMXP = get(gca,'clim');
		end
		xlabel('Frame 1');ylabel('Frame 2');
		DOCUlast = ['2-time correlation averaged over all ROIs'];
		title(char(TITLEstr1,[DOCUlast]));
        prettyplot(gca);
		colorbar
    
    figure;
    set(gcf,'Name','2-time Alternate ave all ROIS');
    HS = pcolor(CA);shading flat
		if ~isempty(CLIMXA);
			set(gca,'clim',CLIMXA);
		else
			CLIMXA = get(gca,'clim');
		end
		xlabel('Frame 1');ylabel('Frame 2');
		DOCUlast = ['Alt 2-time correlation averaged over all ROIs'];
		title(char(TITLEstr1,[DOCUlast]));
			prettyplot(gca);
		colorbar
    
    if 0
    % Average off-diagonal during growth
    if isempty(POINTSUMS); POINTSUMS = [0 length(Xsteps)-ImageJ];end
    Ng = max(POINTSUMS+ImageJ);
    MC = zeros(Ng,1);
    MA = zeros(Ng,1);
    for ii = 1:Ng
        MC(ii) = mean(CC(ii+1:Ng,ii));
        MA(ii) = mean(CA(ii+1:Ng,ii));
    end
    
    figure; % plot of mean off-diagonal vs frame
    set(gcf,'Name','Plot mean off-diagonal to Max POINTSUMS');
    HL = line([1:Ng]-ImageJ,MC);
    xlabel('Frame number');
    ylabel('Off-diag mean corr.');
    title(char(['Plot mean off-diagonal to Max POINTSUMS[',int2str([min(POINTSUMS) max(POINTSUMS)]),']'],TITLEstr1));
    	 	 makeyline((lastframes));
	prettyplot(gca);
    
    figure; % plot of mean off-diagonal vs frame
    set(gcf,'Name','Plot alternate mean off-diagonal to Max POINTSUMS');
    HL = line([1:Ng]-ImageJ,MA);
    xlabel('Frame number'); 
    ylabel('Off-diag mean alt. corr.');
    title(char(['Plot alternate mean off-diagonal to Max POINTSUMS [',int2str([min(POINTSUMS) max(POINTSUMS)]),']'],TITLEstr1));
    	 	 makeyline((lastframes));
        prettyplot(gca);
    end
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(WhatToDo==8); % Look for static speckle after growth  (DO NOT USE YET)
    % sum over POINTSUMS
    %if isempty(POINTSUMS);POINTSUMS = [5 15];end
    if isempty(POINTSUMS); POINTS=[1 length(Xsteps)]-ImageJ;
	else POINTS = POINTSUMS;
	end
	
    Ni = [POINTS(1,1): POINTS(1,2)]+ImageJ;
    Istat = sum(II(:,:,Ni),3);
    Irst = ones(0,1); % get all pixels in ROIs
    for ii = 1:size(ROIS,1);
        Irsti = Istat([ROIS(ii,3):ROIS(ii,4)]+ImageJ,...
				[ROIS(ii,1):ROIS(ii,2)]+ImageJ);
        Irst = [Irst; Irsti(:)];
    end
    
    [Nhist Xhist] = hist(Irst,[0:max(Irst)]);
    
    % Fit to negative binomial distribution
    
    x = Xhist';
    y = Nhist'/length(Irst);
    w = ones(size(x));

    % fit using leasqrs
    global verbose;
    verbose = [0 0];
    stol = 0.000001;
    niter = 100;

    pin = [mean(Irst) 10]; % kbar M
    dp =   [1   1]'*0.0001; % Fraction to vary for deriv

    [ycalc,pout,kvg,iter,corp,covp,covr,stdresid,Z,r2]=leasqr(x,y,pin,'NB_fn',stol,niter,sqrt(w),dp);
    % Run 51 got error in leasqr, Input to SVD must not contain NaN or Inf (called in leasqr)
    
    figure,clf;
    set(gcf,'Name',['Static Speckle Plots']);
		axes('Box','on');
		HL = line(x,y);
		set(HL,'Marker','o','LineStyle','none');
		HLcalc = line(x,ycalc);
    set(HLcalc,'LineStyle','-');
    xlabel('Photons per pixel');
    ylabel('Probability');
    title(char(['Static Speckle Plot : [kbar,M]:[' num2str(pout') ']' ],TITLEstr1));
    prettyplot(gca);
    disp(['kbar, M: ' num2str(pout')])
    	

end
