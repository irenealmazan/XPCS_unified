% everything in here is relevant to creating the sdata (spec data stuff) and the idata (image data stuff)
% The aspects that we might want to replot changed are NOT in here (but note that the way we reduce the image
% data (normalization?, flatfield?, background? would require a new 'redoing' as 

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
	[Iwin,timestampX] = load_MPX3(FullNameHDF,[],badpixfile);
	% The following are not used for tifs, only for pixirad in window mode
	% Note the timestamp for the medipix is only to 1 second, not to milliseconds
	% So need to use the spec information for time.
	% timestamp for medipix also weird in 2017-08
	Icolors = [];
	THRESHuse = [];   %% hardwiring to not use the pixirad high low, etc
	THRESHusedoc = [];
	timestampX_flag = 0;   % timestamp flag from tif is not great, keep using spec
end

		

%------------   If keeping track of valves to plot
%if ~strncmp('caxis',SCNFUNC,min([length('caxis') length(SCNFUNC)]));  % which diffractometer

%	valve_p	= sdata.DATA(:,chan2col(sdata.collabels,'valve_p'));
%	valve_v	= sdata.DATA(:,chan2col(sdata.collabels,'valve_v'));
	valve = sdata.DATA(:,chan2col(sdata.collabels,valvename));
	lastframes = find(diff(valve)>=diffvalve);  % this is automatically in spec points, when using to index matlab point +1
	disp(['Last frames (Spec point number convention) at valve switches: ' num2str(lastframes')]);
 
		% we end up always zeroing to the first point in scan, but if we need to compare between scans
		if timestampX_flag==0;
			%timestampX = timestampSpec + timestampEpoch(1);
			timestampX = timestampSpec;
		end
	



	

%	Here we are being consistent with enumeration in SPEC points and ImageJ
%		ImageJ convention is used in the SPEC ROI's and EPICS ROIs images
%		SPEC convention is used in its points when we 'write' down a data point
%   We found various indices in matlab, if ImageJ=1, then the following will turn them
%		into spec/imageJ indexing, if ImageJ=0, then we remain using matlab indexing
	SPECpts = [1:(length(timestampX))]		- ImageJ;
	YROWpts = [1:(length(Iwin(:,1,1)))]	- ImageJ;  
	XCOLpts = [1:(length(Iwin(1,:,1)))]	- ImageJ;
	
	
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

% We are not using the Iwinnormb (from the pixirad? in windowed version which we haven't used in long time so I did not keep up
