% called by datainput_AREA which set ups various initializations, and also uses a sample file
% took out fixlogimage since matlab past 2016? will not croak in an image or pcolor
% 	if one of the pixels is NaN (or 0) (and it seems to do it sensibly by ignoring them)

clear Iwin Iwinnorm;
HLine=[];HSurf=[];

if HDF5==1   % not really used, but legacy for when using pixirad with windows
	STR.scanflag=1;
	STR.imname=imname;
	STR.SPECpath=SPECpath;
	STR.IMAGEpath=IMAGEpath;
	STR.p_image=imagepointnum;
	STR.ending='.h5';
else  % medipix, pilatus, all saved in tif thus far (usually)
	STR.scanflag=1;
%	STR.imname=[];      % image names use standard unique naming schemes
%	STR.imname='pil_';  % caxis non-unique naming scheme for images (sigh) 
% 	STR.p_image = 'p_image';
	STR.SPECpath=SPECpath;
	STR.IMAGEpath=IMAGEpath;   % need to change pathdisplay if need pilatus
	STR.p_image=imagepointnum;   %% STR.p_image=[] if it uses the preferred points instead
	STR.imname = imname;
	STR.ending='.tif';
end

% If ImageJ = 1; (then we are using ImageJ and spec indexing in ROI's and scan points (start at 0)
%			= 0; Then we are using matlab indexing in ROI's and scan points
%				Note - most functions take the indexing as given (or will have an ImageJ flag)
%				Those that take the indexing as given should be given NDX + ImageJ and they will
%				work no matter which (since ImageJ will take account of it, 0 (no change) +1 change)
ImageJ = 1; % true - using imageJ/spec indexing [start at 0] 

% FORCE only one scan at a time
SCNs = eval(SCNstr); 
SCNs = SCNs(1); 

%if strncmp('caxis',SCNFUNC,min([length('caxis') length(SCNFUNC)]));  % which diffractometer
%	INFOinterest = ['Sam_Z  wcay  hrmy  stemp  TABLEVERT  Delta  Chi  Eta  s6vgap  s6hgap'];     %!!! put 2 spaces between each!!!!
%else
%	INFOinterest = ['s6vgap  s6hgap  Phi  TZ  TY  Chi  Eta  Delta  Nu  TX'];     %!!! put 2 spaces between each!!!!
%end

%----------  Make the matrix of the names for the images, note sdata holds information from spec file
[NameMatrix,sdata]	= make_imnames_2018_08(HDFfilename,SCNs,STR);  % added multiple file support for 2018_08
FullNameHDF 		= addnames2matrix([STR.IMAGEpath,filesep],NameMatrix.fullfilenames);

%----------  At this point can pull out quite a bit from the spec file (in sdata)

	% with make_imnames_2018_03 which permits more than one scan to be concatenated
	%	the scnheaers and the multiple complex stringgs send up being in different cells
	%	for now, assume no changes in any other positions that we use here
	%	but in future, should incorporate this and permit between scans, etc
	%  
	if iscell(sdata.scnheader);    % cat won't do anything much with strings, if they are different so I'm not quite
		sdata.scnheader		= cat(1,sdata.scnheader{1});
		sdata.fileheader	= cat(1,sdata.fileheader{1});
		sdata.comments		= cat(1,sdata.comments{1});
		sdata.SCNTYPE		= cat(1,sdata.SCNTYPE{1});
		sdata.SCNDATE		= cat(1,sdata.SCNDATE{1});
	end
	
	% I expect that I will make readspechelp files compatible with cells for the above
	[OM,scninfo,OMspecial]	= readspechelp(sdata.scnheader,sdata.fileheader,sdata.comments);
	INFOndx 				= chan2col(scninfo.allangles_label,INFOinterest);
	INFOpositions 			= scninfo.allangles_i(INFOndx);
		DUM=[];
		% put into a line
		for ii=1:length(INFOpositions);
			DUM = [DUM,num2str_round(INFOpositions(ii)) ' : '];
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
	meansecperpoint = mean(secperpoint);

	%---- get the normalization corrections ready
	% this set was designed for more or less constant time but when 
	% slits were changing a lot (to make intensities more or less on
	% similar scales for purposes of plotting
	
	%---- get the normalization corrections ready, these are designed 
	%---- when counting to monitor, (time can fluctuate wildly in some cases)
	%---- This gives something scaled back to 'raw' counts
	if NORMFLAG == 1;	% using = 1 (det*filters/mon ,*monave)  = 2 (det/mon)  = 3 (raw)
		Norm = MONave .*  filter_correction(filterpoint,filter1,filtercorr)./(Monitor(:,1));
		NORMdoc = ['Norm = monave*filtercorrection/(Monitor per point): WITH filter correction'];
	elseif NORMFLAG ==2; 
		Norm = MONave ./(Monitor(:,1));
		NORMdoc = ['Norm = monave/(Monitor per point): No filter correction'];
	elseif NORMFLAG ==3;
		Norm = ones(Monitor(:,1));
		NORMdoc = ['Norm = no corrections'];
	else
		Norm = MONave ./(Monitor(:,1));	
		NORMdoc = ['Norm = monave/(Monitor per point): No filter correction'];
	end
	
	disp(char('* *','Normalization that is used is as follows',NORMdoc,'* *'));


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
		[pfilename(sdata.SCNTYPE)],...
		[pfilename(DOCU0),' ',pfilename(DOCUscan)]);

TITLEstrshort = char(...
		[pfilename(HDFfilename),' #',SCNstr,' : ',DOCUInt,' [',pfilename(DOCU0),']'],...
		[pfilename(sdata.SCNTYPE)]);
		

%------------   If keeping track of valves to plot
%if ~strncmp('caxis',SCNFUNC,min([length('caxis') length(SCNFUNC)]));  % which diffractometer

%	valve_p	= sdata.DATA(:,chan2col(sdata.collabels,'valve_p'));
%	valve_v	= sdata.DATA(:,chan2col(sdata.collabels,'valve_v'));
	valve = sdata.DATA(:,chan2col(sdata.collabels,valvename));
	lastframes = find(diff(valve)>=diffvalve);  % this is automatically in spec points, when using to index matlab point +1
	disp(['Last frames (Spec point number convention) at valve switches: ' num2str(lastframes')]);
%else valve = sdata.DATA(:,chan2col(sdata.collabels,valvename));
%	lastframes = find(abs(diff(valve))>=diffvalve);
%	disp(['Last frames (Spec point number convention) at valve switches: ' num2str(lastframes')]);
%end  
	
		if timestampX_flag==0;
			timestampX = timestampSpec + timestampEpoch(1);
		end
	

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
	
	if strncmp(SCNXLABEL,etimename,4); %usually the elapsed time name is 'Time'
		SCNXLABEL	= [SCNXLABEL,' [\Delta sec] from SpecFile '];  % in pixirad, 
		Xsteps	= timestampX-timestampX(1);
	elseif Xstepcol<1;
		Xsteps	= [1:length(timestampX)]-ImageJ;     
	else
		Xsteps	= sdata.DATA(:,Xstepcol);
	end

	

%	Here we are being consistent with enumeration in SPEC points and ImageJ
%		ImageJ convention is used in the SPEC ROI's and EPICS ROIs images
%		SPEC convention is used in its points when we 'write' down a data point
%   We found various indices in matlab, if ImageJ=1, then the following will turn them
%		into spec/imageJ indexing, if ImageJ=0, then we remain using matlab indexing
	SPECpts = [1:(length(Xsteps))]		- ImageJ;
	YROWpts = [1:(length(Iwin(:,1,1)))]	- ImageJ;  
	XCOLpts = [1:(length(Iwin(1,:,1)))]	- ImageJ;
	
	Nr = length(YROWpts);
	Nc = length(XCOLpts);
	Nz = length(SPECpts);
	
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
			%disp('BKG empty');
			IIimbkg = 0;
			Iwinimbkg = 0;
		elseif length(BKG)>1
			IIimbkg		= mean(IInorm(:,:,[BKG(1):BKG(2)]+ImageJ),3);
			Iwinimbkg	= mean(Iwinnorm(:,:,[BKG(1):BKG(2)]+ImageJ),3);
		else  % test of this other detector background stuff
			BKGshoulders = [20 60;160 170];  % 
			%BKGshoulders = [20 80;140 170];  % 
			if ~isempty(lastframes)
			%IIimbkg=0;
		    IIimbkg =  image_est_bkg(IInorm(:,:,[2:lastframes(1)]+1),'Y',BKGshoulders);  % pre 2018_0315_1
		    else
		    %IIimbkg=0;
		    IIimbkg =  image_est_bkg(IInorm,'Y',BKGshoulders);
		    end
		    Iwinimbkg = 0;
		    DOCUInt = [DOCUInt(1:end-1) '-B)'];
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
%		HS = pcolor(XCOLpts,YROWpts,(IInormb(:,:,SINGLE+ImageJ)));
		HS = imagesc(XCOLpts,YROWpts,(IInormb(:,:,SINGLE+ImageJ)));
		shading interp
		HSurf.IInormb = HS;
		
			
			if ~ASPECTflag
			axis equal;
			else
			HA = gca;
			daspect(gca,ASPECT);   %this might not be in earlier version
			HA = gca;,
			set(HA,'DataAspectRatioMode','manual','DataAspectRatio',[ASPECT]);
			end
			
			if ~isempty(AXISdet); axis(AXISdet);end
			if ~isempty(CLIMlin); 
				set(gca,'clim',CLIMlin);
				DOCUclim = mDOCUclim(CLIMlin);
			end;
	else
		figure;clf
		set(gcf,'Name',['Log10 Image',THRESHusedoc,' Spec #' int2str(SINGLE)]);  % assumes SINGLE in spec/ImageJ
%		HS = pcolor(XCOLpts,YROWpts,log10(fixlogimage((IInormb(:,:,SINGLE+ImageJ)))));
		HS = imagesc(XCOLpts,YROWpts,log10((IInormb(:,:,SINGLE+ImageJ))));
		shading flat
		HSurf.IInormb = HS;

			
			if ~ASPECTflag
			axis equal;
			else
			HA = gca;
			daspect(gca,ASPECT);   %this might not be in earlier version
			HA = gca;,
			set(HA,'DataAspectRatioMode','manual','DataAspectRatio',[ASPECT]);
			end
			
			
			if ~isempty(AXISdet), axis(AXISdet);end
			if ~isempty(CLIMlog); 
				set(gca,'clim',CLIMlog);
				DOCUclim = mDOCUclim(CLIMlog);
			end;
	end
	
	DOCUlast = ['(Spec pt #',int2str(SINGLE),') ', DOCUInt, DOCUclim];
	title(char(TITLEstr1,DOCUlast,pfilename(INFOstr)));
	xlabel(XCOL);ylabel(YROW);
	%prettyplot(gca,[2 2 4 3.8]);
	plot_adjust(gca,[0 1],[],[2 2 4 3.8]);
	colormap(COLORMAP);
end	


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(WhatToDo==1) %& ~isempty(ROIS) %% open up make ROIS
	myROIS = makerois(4,gcf);
	xlabel(XCOL);ylabel(YROW);  % must have done the WhattoDo =0
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(WhatToDo==2)

	[HROI,COLORORDER] = showrois(ROIS,gcf);
		
		figure;clf;
		set(gcf,'Name','Enumeration of ROIs colors');

%		prettyplot(gca);    % put here or it clobbers the showROIS colororder 		
		showrois(ROIS,gcf,2,COLORORDER);
		title(char('Enumeration of the ROIs and their colors used for',...
			[pfilename(HDFfilename),' #', SCNstr,' : ', sdata.SCNDATE]));

		axis equal;
		if ~isempty(AXISdet); 
			axis(AXISdet);
		else
			axis([min(XCOLpts) max(XCOLpts) min(YROWpts) max(YROWpts)]);
		end

			LEGEND = [char([ones(length(HROI),1)*'ROI #']) int2str([1:length(HROI)]')];
			legend(LEGEND);
			xlabel(XCOL);ylabel(YROW);  % must have done the WhattoDo =0
			
			plot_adjust(gca,[0 1]);
		
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
	if ~LOGFLAG
		HS = imagesc(XCOLpts,YROWpts,((IInormb(:,:,SINGLE+ImageJ))));  
	else
		HS = imagesc(XCOLpts,YROWpts,log10(IInormb(:,:,SINGLE+ImageJ)));  
	end
	
	colormap(COLORMAP)
	shading flat
		HT = title(...
			char([TITLEstrshort],...
			['Spec Pt # [ 00000 ] : step =  [0000000]',mDOCUclim(CLIMlog)]));	

		if ~ASPECTflag
		axis equal;
		else
			HA = gca;
			daspect(gca,ASPECT);   %this might not be in earlier version
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
		
		if ~isempty(CLIMlog) & LOGFLAG; 
			set(gca,'clim',CLIMlog);
		elseif ~isempty(CLIM) & ~LOGFLAG;
			set(gca,'clim',CLIMF);
		end;
		plot_adjust(gca,[0 1]);
		
		if any(WhatToDo==2);
			if exist('COLORORDER','var')%length(ROIS(:,1))<length(COLORORDER(:,1))
				showrois(ROIS,gcf,0.5,COLORORDER);
			else
				[HROI,COLORORDER]=showrois(ROIS,gcf,0.5);
			end
		end		
		

	% Prepare the new file.
		if exist('M1.avi')==2; delete('M1.avi'); end  
		vidObj = VideoWriter('M1.avi');
		vidObj.FrameRate = 5;		% default 30 fps


		open(vidObj);


%	for ii=1:length(timestampX)  % CT
	for ii=1:length(Xsteps);
		PTstr = int2str(ii-1+10000);PTstr = PTstr(2:end);
		
		if ~LOGFLAG
			set(HS,'CData',((IInormb(:,:,ii))));
			if ~isempty(CLIMlin); set(gca,'clim',CLIMlin);end;
		else
			set(HS,'CData',(log10(IInormb(:,:,ii))));
			if ~isempty(CLIMlog); set(gca,'clim',CLIMlog);end;
		end		
		set(HT,'String',...
			char([TITLEstrshort],...
			['Spec Pt # [',PTstr,'] : step = [',num2str(Xsteps(ii)),']',mDOCUclim(gca)]));
		currFrame = getframe(gcf);
		writeVideo(vidObj,currFrame);
	end
	
	% have short pause at end of movie while sitting on the end point
	for ii=1:3;
		currFrame = getframe(gcf);
		writeVideo(vidObj,currFrame);
	end
	
	disp('*');
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

% N_ROIS_SUM vector with the ROIS to use, SoX(Y)FLAG whether to sum in X or Y for that ROI
if length(SoXFLAG)~=length(N_ROIS_SUM);
	SoXFLAG = ones(1,length(N_ROIS_SUM));
end
if length(SoYFLAG)~=length(N_ROIS_SUM);
	SoYFLAG = ones(1,length(N_ROIS_SUM));
end	

iii = 0;
for ii=N_ROIS_SUM; 
  iii=iii+1;
		
  
  if SoXFLAG(iii)
	figure;clf
	set(gcf,'Name',['SoX ROIs #' int2str(ii)]);

	 if ~LOGFLAG;
		HS = imagesc(Xsteps,SoX.ndx{ii},SoX.images{ii});
%		HS = surf(Xsteps,SoX.ndx{ii},SoX.images{ii}, SoX.images{ii}));
	 else
		HS = imagesc(Xsteps,SoX.ndx{ii},log10(SoX.images{ii}));
%		HS = surf(Xsteps,SoX.ndx{ii},log10(SoX.images{ii}),log10(SoY.images{ii}));
	 end

	 HPC.SoX_ROIS(ii) = HS;
		shading flat; view(0,90); colormap(COLORMAP)
		xlabel(SCNXLABEL);ylabel(YROW);
		title(char(TITLEstr1,[DOCUInt, ' summed over XCOL in ROI #' int2str(ii)]));
	 set(gca,'Xlim',[min(Xsteps) max(Xsteps)]);
	 
		% Now YLim is simply Y on the graph
		%	if ~isempty(AXISdet); set(gca,'Ylim',AXISdet(3:4));end
		if diff(ROIS(ii,[3:4]))   % if non zero do this
			set(gca,'Ylim', ROIS(ii,[3:4])); 
	 	else % accommodate line ROI that user summs in non-useful dirction
			set(gca,'Ylim', ROIS(ii,[3])+[-5 5]);
		end

	 			
	plot_adjust(gca,[0 1]);
	 
		
		if ~isempty(lastframes)
			HYline = makeyline(Xsteps(lastframes+ImageJ),'b',gca);  % and makes a straight line
		end

  end
	
  if SoYFLAG(iii)  % don't plot if SoY flag 0
	figure;clf
	set(gcf,'Name',['SoY ROIs #' int2str(ii)]');
	
	 if ~LOGFLAG;  % surf(X,Y,Z,Color)
%	 HS = surf(Xsteps,SoY.ndx{ii},((SoY.images{ii})),((SoY.images{ii})));	
	 HS = imagesc(Xsteps,SoY.ndx{ii},((SoY.images{ii})));		
	 else
%	 HS = surf(Xsteps,SoY.ndx{ii},log10(fixlogimage(SoY.images{ii})),log10(fixlogimage(SoY.images{ii})));
	 HS = imagesc(Xsteps,SoY.ndx{ii},log10(SoY.images{ii}));
	 end

	 HPC.SoY_ROIS(ii) = HS;
	 shading flat; view(0,90); colormap(COLORMAP)
	 xlabel(SCNXLABEL);ylabel(XCOL);
	 title(char(TITLEstr1,[DOCUInt ' summed over YROW in ROI #' int2str(ii)]));
	 set(gca,'Xlim',[min(Xsteps) max(Xsteps)]);

		% Now YLim is simply Y on the graph
		%	if ~isempty(AXISdet); set(gca,'Ylim',AXISdet(3:4));end
		if diff(ROIS(ii,[1:2]))   % if non zero do this
			set(gca,'Ylim', ROIS(ii,[1:2])); 
	 	else % accommodate line ROI that user summs in non-useful dirction
			set(gca,'Ylim', ROIS(ii,[1])+[-5 5]);
		end



	plot_adjust(gca,[0 1]);

		if ~isempty(lastframes)
			HYline = makeyline(Xsteps(lastframes+1),'b',gca);  % and makes a straight line
		end

  end
end

	if  1    %   plots that show line plots as function of nu or del
		% Note - if any SoX (SoY) requested  on any ROIs, the all ROIS will done
		if any(SoXFLAG)
		figure; clf
		set(gcf,'Name','SoX and Sum Points ');
		HL = semilogy(SoX.ndx{1},sum(SoX.images{1}'));
					% in event only one point because of dimensionality
					if size(sum(SoX.images{1}'))<2;HL(1).Marker='o';end
		
			for jj = 2:length(SoX.images);
				HL(jj) = line(SoX.ndx{jj},sum(SoX.images{jj}'));
					% in event only one point because of dimensionality
					if size(sum(SoX.images{jj}'))<2;HL(jj).Marker='o';end
			end

		xlabel(YROW),ylabel('Int (arb)'),
		title(char(TITLEstr1,['summed over XCOLS in all ROI and ',int2str(SIZE3),' scan pts']));
		plot_adjust(gca,0); % 
			for ii=1:length(ROIS(:,1));
				set(HL(ii),'Color',COLORORDER(ii,:),'LineWidth',2,'DisplayName',['ROI #',int2str(ii)]);
			end
				Hlegend = legend('show','AutoUpdate','off');

		
		
		HLine.SoX_SPts = HL;clear HL;
	
		end
	
		if any(SoYFLAG)	
		figure; clf
		set(gcf,'Name','SoY and Sum Points ');
		HL = semilogy(SoY.ndx{1},sum(SoY.images{1}'));
							% in event only one point because of dimensionality
					if size(sum(SoY.images{1}'))<2;HL(1).Marker='o';end
			for jj = 2:length(SoY.images);
				HL(jj) = line(SoY.ndx{jj},sum(SoY.images{jj}'));
					% in event only one point because of dimensionality
					if size(sum(SoY.images{jj}'))<2;HL(jj).Marker='o';end
			end

		xlabel(XCOL),ylabel('Int (arb)'),
		title(char(TITLEstr1,['summed over YROWS in all ROI and ',int2str(SIZE3),' scan pts']));
		plot_adjust(gca,0); %
			for ii=1:length(ROIS(:,1));
				set(HL(ii),'Color',COLORORDER(ii,:),'LineWidth',2,'DisplayName',['ROI #',int2str(ii)]);
			end
				Hlegend = legend('show','AutoUpdate','off');
	 
		HLine.SoY_SPts = HL;clear HL;
		end


		if ~isempty(POINTSUMS)   % sum over points sets in POINTSUMS, one plot per ROI
	
		jjj=0; 
	
			for jj=N_ROIS_SUM;
				jjj=jjj+1;
				
				idROI = int2str(jj);
				if SoXFLAG(jjj)
		
				figure; clf;
				set(gcf,'Name',['SoX ROI # ' int2str(jj) ' and Sum select points']);
					for ii=1:length(POINTSUMS(:,1));
					Ni = [POINTSUMS(ii,1): POINTSUMS(ii,2)]+ImageJ;  %use as matlab
%					HL(ii) = line(SoX.ndx{1},sum(SoX.images{1}(:,Ni)'));
					HL(ii) = line(SoX.ndx{jj},sum(SoX.images{jj}(:,Ni)'./length(Ni)));
					end
		
		
				set(gca,'Yscale','log')
				xlabel(YROW),ylabel('Int (arb)'),
				title(char(TITLEstr1,['# ' int2str(jj) ' ROI summed over X, then summed over Spec Point Range']));
				plot_adjust(gca,0); %
					legend(addnames2matrix('sumX between [', int2str(POINTSUMS),']/Npts'))
				HLine.SoX_SoSelectPt = HL;clear HL;
				set(HLine.SoX_SoSelectPt,'linewidth',1.5);
				end
		
				if SoYFLAG(jjj)

				figure; clf;
				set(gcf,'Name',['SoY ROI # ' int2str(jj) ' and Sum selected points']);
					for ii=1:length(POINTSUMS(:,1));
						Ni = [POINTSUMS(ii,1): POINTSUMS(ii,2)]+ImageJ;
%						HL(ii) = line(SoY.ndx{1},sum(SoY.images{1}(:,Ni)'));
						HL(ii) = line(SoY.ndx{jj},sum(SoY.images{jj}(:,Ni)'./length(Ni)));
					end
				set(gca,'Yscale','log')
				xlabel(XCOL),ylabel('Int (arb)'),
				title(char(TITLEstr1,['# ' int2str(jj) ' ROI summed over Y, then summed over Spec Point Range']));
				plot_adjust(gca,0); %
					legend(addnames2matrix('sumY between [', int2str(POINTSUMS),']/Npts'));
				HLine.SoY_SoSelectPt = HL;clear HL;
				set(HLine.SoY_SoSelectPt,'linewidth',1.5);
				end
			end   % end of if jj % jj=N_ROIS_SUM;
	
		end  % sum over points sets in POINTSUMS
	end %%   plots that show line plots as function of nu or del
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% note if the ROI has NaN in any of the pixels, this will mess up and give NaN
%%% need to fix sumrois to fix that
if any(WhatToDo==5)  %% make summed lines of ROI's on images
	figure;clf
	set(gcf,'Name','ROIs as counters');
	 
	 [sumROIS,sumROISnormed] = sumrois(IInormb,ROIS,ImageJ);
	 %  sumROIS sums over the ROI, sumROISnormed divides by the number of pixels
%	 [sumROIS,sumROISnormed] = sumrois(IInormb,ROIS,ImageJ);
	 %	Isum	Array  : each column for each ROI. Length of vectors (rows) num of images
	 %HL = semilogy(Xsteps,sumROIS);
	 HL = plot(Xsteps,sumROISnormed); 
		YLABEL = '[Inorm(summed over ROI)/ROIsize]';
%	 HL = plot(Xsteps,log10(sumROISnormed)); YLABEL = 'log10[Inorm(summed over ROI)/ROIsize]';


	 xlabel(SCNXLABEL);ylabel(YLABEL);
	 title(char(TITLEstr1,pfilename(INFOstr),YLABEL));
	 plot_adjust(gca,0); %
% mocvd
%	 makeyline(Xsteps(lastframes+ImageJ),'b',gca);   % want colored lines on this plot, not white	 
	 for ii=1:length(ROIS(:,1));
		set(HL(ii),'Color',COLORORDER(ii,:),'LineWidth',2,'DisplayName',['ROI #',int2str(ii)]);
	 end

	 Hlegend = legend({HL.DisplayName},'AutoUpdate','off');   % why does it work here and not in 4?
	 makeyline(Xsteps(lastframes+ImageJ),'b',gca);   % want colored lines on this plot, not white	 
	 
	 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(WhatToDo==6);   %to make it automatic for summs   % ? is Xsteps missing one point
	if isempty(POINTSUMS); 
		POINTS=[1 : length(Xsteps)]-ImageJ;
	else 
		POINTS = POINTSUMS;
	end

for ii = 1:length(POINTS(:,1))
	Ni = [POINTS(ii,1): POINTS(ii,end)]+ImageJ;
	NP = length(Ni);
	POINTstr = ['[ ' int2str(Ni(1)-ImageJ) ':' int2str(Ni(end)-ImageJ) ']'];
	
	figure;clf;
	set(gcf,'Name','Image summed over points/Num Points');	
	if ~LOGFLAG

		HS = imagesc(XCOLpts,YROWpts,sum(Iwinnorm(:,:,Ni),3)./(NP));
		colormap(COLORMAP)
			if ~ASPECTflag
				axis equal;
			else
				HA = gca;
				daspect(gca,ASPECT);   %this might not be in earlier version
				HA = gca;,
				set(HA,'DataAspectRatioMode','manual','DataAspectRatio',[ASPECT]);
			end

			if ~isempty(AXISdet); axis(AXISdet);end		
			if ~isempty(CLIMlin); 
				set(gca,'clim',CLIMlin);
				DOCUclim = mDOCUclim(CLIMlin);
			else 
				CLIMtemp = get(gca,'clim');
				DOCUclim = mDOCUclim(CLIMtemp);
			end;
			DOCUInt6 = [DOCUInt, ' SUMMED over SPEC scan pts ',POINTstr];
	else
	
		HS = imagesc(XCOLpts,YROWpts,log10(sum(IInormb(:,:,Ni),3)./(NP)));
		colormap(COLORMAP)	
			if ~ASPECTflag
				axis equal;
			else
				HA = gca;
				daspect(gca,ASPECT);   %this might not be in earlier version
				HA = gca;,
				set(HA,'DataAspectRatioMode','manual','DataAspectRatio',[ASPECT]);
			end
			
			if ~isempty(AXISdet), axis(AXISdet);end
			if ~isempty(CLIMlog); 
				set(gca,'clim',CLIMlog);
				DOCUclim = mDOCUclim(CLIMlog);
			else 
				CLIMtemp = get(gca,'clim');
				DOCUclim = mDOCUclim(CLIMtemp);
			end;
			DOCUInt6 = [DOCUInt, ' SUMMED over SPEC scan pts ', POINTstr];
	end
	
	shading flat;
	
	if FLATFIELDFLAG
		DOCUInt = [DOCUInt,' [FF]'];
	end
	
	DOCUlast = char(DOCUInt6,pfilename(INFOstr));
	title(char(TITLEstr1,DOCUlast));xlabel(XCOL);ylabel(YROW);
	plot_adjust(gca,0,[],[2 2 4 3.8]); %
		% only put ROIs on if we did something like option 2
		if any(WhatToDo==2);
			if exist('COLORORDER','var')%length(ROIS(:,1))<length(COLORORDER(:,1))
				showrois(ROIS,gcf,1.5,COLORORDER);
			else
				[HROI,COLORORDER]=showrois(ROIS,gcf,1.5);
			end
		end
	
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(WhatToDo==7);  % plot flatfield FLATFIELDFLAG

	figure
	set(gcf,'Name',['Flat Field']);
	HS = imagesc(XCOLpts,YROWpts,FLATFIELD);
	shading flat
	axis equal;
	colormap(COLORMAP)
    if ~isempty(AXISdet); axis(AXISdet);end
	title(['Flat Field: from ' pfilename(FlatField.filenameroot) ' [' FLATFIELDinfo ']']);
	xlabel(XCOL);ylabel(YROW);
	plot_adjust(gca,0,[],[2 2 4 3.8]); %
    %set(gca,'position',[.17 .1 .65 .68]);

			if exist('COLORORDER','var')%length(ROIS(:,1))<length(COLORORDER(:,1))
				showrois(ROIS,gcf,1,COLORORDER);
			else
				[HROI,COLORORDER]=showrois(ROIS,gcf,1);
			end	
 
end	

%---------------------  GBS sections WhattoDo 8 and above
 % user_input_guts_2timesection
%---------------------------




