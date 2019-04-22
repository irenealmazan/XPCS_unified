function [idata,sdata] = datainput_AREA(user_sample_AREA,ReUseFlag)
% for now, just a calling function to keep the workspace cleaner
% Uses user_input (the 'user thing to edit  (>> datainput_AREA('samplefile')
% which itself calls user_input_guts_XXX  which is a big horrendous program
% CT 2016-10 changed to POINTSSUMS and SINGLE refer to points in SPEC
%		convention instead of matlab (spec starts at 0)
%		earlier in 2016-10 it would convert in guts
%		and even earlier, here and in guts used matlab point number

SoX=[];SoY=[];  
RUNPARAMS=parameterfiledisplay('showfilename');

badpixfile = 'badpix_CTpil_2018_08.txt';   % use this bad pix file and replace bad pix with NaN
%badpixfile = [];  % don't replace bad pixels with NaN

if nargin<1;
	user_sample_AREA = 'user_datainputAREA';
end	
if nargin < 2;
	ReUseFlag = 0;
end


%% Do NOT touch these, they set up empties used for defaults
[SPECpath,IMAGEpath,COMMONpath,HOMEpath] = pathdisplay;
NORMFLAG=0;
% From FORCE.FLAGdet		= 1 (det*filters/mon)  = 2 (det/mon)  = 3 (raw use default - which is /sec*mon then scaled back to average mon
ASPECTflag=0;   % see daspect if ASPECTflag=1; also need ASPECT
WRITETIFflag = 0;
COLORMAPflag = 0;
COLORMAP = 'parula';  % matlab default   (use as colormap(COLORMAP)
COLORMAP = 'jet';

%  example of some extra information about state of motors at beginning of run brought out make sure to put 2 spaces between each
%INFOinterest = ['Sam_Z  wcay  hrmy  stemp  TABLEVERT  Delta  Chi  Eta  s6vgap  s6hgap'];     %!!! put 2 spaces between each!!!!
%INFOinterest = ['s6vgap  s6hgap  Phi  TZ  TY  Chi  Eta  Delta  Nu  TX'];     %!!! put 2 spaces between each!!!!
INFOinterest = [];

BKG=[];  %=1 for the quick background before valve fit and use
		% = [num1 num2] to use average of image between frames num1 and num2

CLIM=[];
CLIMlin=[];
CLIMlog=[];
AXISdet = [];   %Xmin Xmax Ymin Ymax
ROIS=[];
POINTSUMS = [];
SINGLE=0;  % (if no other choice, plot the 1st image, which is pt zero in spec) 
THRESHuse=[];
FLATFIELDFLAG = 0;  % 1 to use a flatfield (need to specify a mat file where it exists)
detwallFLAG = 0;  % used for different rotations of detector on wall and on arm

POINTSUMS = [];
SoXFLAG = 1;
SoYFLAG = 1;
N_ROIS_SUM = [1];    % all ROIS are summed for the ROI as 'detector', but which ones
					% are summed over Y (or over X) for the plots of space over scan
					% default was the 1st ROI, this allows several or a diffeent one


fixfile(RUNPARAMS);	%during run this is usually neccesary to re-index updated spec files, 
					% RUNPARAMS has an enddate so this will ding if past date

eval(RUNPARAMS);  % note this will set a variety of things - they can be overridden by adding them 

%	in the sample file called below
%% Go into the file as named below for user changes (poor mans gui)
eval(user_sample_AREA);  %<< can change to a different file, but see this file for entry


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%   For caxis, moved some of the stuff in the datainput user file for a sample here
%%%%%%%   since it was likely not to change, and to make it easier to 
%%%%%%%   have anyone edit the user file per sample and get OK results
COLORMAPflag = 1;
%FORCE.FLAGdet = 2;  % = 1 (det*filters/mon)  = 2 (det/mon)  = 3 (raw 

%AXISdet = [1 487 1 195]-1;   % to put it into spec like

MOVIEROIFLAG = 0;  % [1] use ROI#1 to frame movie, [0] use AXISdet (if defined) or whole image'

ASPECTflag=0;  % do not do any forshortening based 
%ASPECTflag=1;ASPECT = round([1./sind(2.42) 1 1]); % use a forshortening based on q
%ASPECTflag=1;ASPECT = round([1 1./sind(2.4) 1]); % use a forshortening based on q back wall
%			 to calculate on specular type rods, use 1./sind(2theta/2) on the X or Y that is delta 
%%%%%%%  ROIS and AXISdet (using ImageJ (SPEC) pixel conventions  %%%%%
%%%%  that is, start from [0,0] indexing
%%%% Medipix [516x516 pixels] 0 515 1 515   or [1 516 1 516]-1   
%%%% Pilatus 0 486 0 194   or [1 487 1 195]-1

%SoXFLAG = 1;  %% =1 will sums over delta in detector
%SoYFLAG = 1;  %% = 1 will sum over nu in Pilatus detector

%% if you desire to have the 'summed' plots only over a particular set of points
%POINTSUMS = [0 10;SINGLE-5 SINGLE+5]; % scan 22
%POINTSUMS=[];

					%% optional for  [4 5 6 8]
					%% POINTSUMS are [i1 j1;i2 j2;...] where each row are
					
%% IN XPCS, tend to go back and forth, for others, this can typically be set and forgot
%% detwallFLAG = 0;  % detector on detector arm (X is del and Y is nu)  (also typical Pilatus
%% detwallFLAG = 1;  % detector on wall (X is nu and Y is del)
detwallFLAG = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


FLATFIELDFLAG = 0;  %=1 then try applying flatfield correction to images and data constructed from them
% These are related to the 2-time stuff		
CLIMXP = [];
CLIMXA =[];
%---------- no longer needed in the user_sample_AREA choices

		%% applied during [-3 -2 -1]
%Nharm 		= 8;	%% << hand estimated for harmonic contamination  
					%%			in Win - High/Nharm	

		%% applied during [0 3 4 5 6 8]										   
%BKG			= [];	%% =[] no background subtracted,
					%% = [n1 n3] mean images n1:n3 (if using ImageJ indexing for ROI, etc), (1st scan is n1=0)
%BKG = 1;			%%  = do this bit of background from stuff prior valve open					
THRESHuse 	= [];	%% takes value of 'low';'high';'win' ('win' default)
%					%%  
%---------------------------------------

%%%%%%%%%%%%% Detector on Arm or on Back wall? %%detwallFLAG=1 on back wall%%%%%%%%%%%%%%%%%%
if detwallFLAG;
XCOL = 'pix(nu)';YROW = 'pix(del)';DOCUclim=[];  % medipix on back wall rotated, 
else
XCOL = 'pix(del)';YROW = 'pix(nu)';DOCUclim=[];   % medipix on det arm, also, also Pilatus typically mounting
end

%%%%%%%%%%%%%%%%%%%  BELOW RARELY CHANGES BUT COULD %%%%%%%%%%%%%%%%%%%%%%%
%LOGFLAG = 0;   % [1] plot Intensity log ; [0] linear
%NORMFLAG = 1;		% [1] I=I/hexmon; [0] no normalization % I think hardwired so doesn't use this

%%%%%%%%% BELOW CHANGES EVENMORE RARELY BUT MAY NEED TO CHANGE AT TIMES	
hex100mA = 1; %300000;   % 3um x 3um ?, 10um x 3um ?  %% changed - later it will go to the runparams
% for normalization, this should change if s6gap is changed  but we are keeping same
LOWHIGH =   [1 2];   % who is [low high] lo threshold is 1st image, high threshold is 2nd image
HDF5 = [0];	% [1] hdf file; [0] tif file  % haven't done much with this requires some other changes


%%%%%%%%%%%%%% DO NOT TOUCH BELOW %%%%%%%%%%%%%%%%%%%%
% if flatfield exists (use make_flatfield and use the field flatfieldscaled

if FLATFIELDFLAG==1;
	load FlatField_2017_0809_2_num64;
	FLATFIELD = FlatField.imnormnan;
	FLATFIELDinfo = 'imnormnan';
	Badimage = FlatField.Badimage;
	
	% Fix up extra bad pixels  (these are the 'lines')
	% Note - at present, these are not getting back into the Flatfield structure
	% eventually want to fix this
	ibpe = [261 261]; % y values
	jbpe = [176 250]; % x values
	ibpe = [251 261 261 261 277 298 372 388 409 494 513 256   3 261 256]; % y values
	jbpe = [  7  55 140 240 256 256 256 256 256 256 256 280 393 424 499]; % x values

	for ii = 1:length(ibpe)
		Badimage(ibpe(ii),jbpe(ii)) = 1;
		imnormnan(ibpe(ii),jbpe(ii)) = NaN;
		FLATFIELD(ibpe(ii),jbpe(ii)) = NaN; % GBS added this, above line did not seem right
	end
	
else
	FLATFIELD=1;
	BadImage=0;
end

if LOGFLAG
	CLIMlog = CLIM;
	CLIMlin = 10.^(CLIMlog);
else
	CLIMlin = CLIM;
	if ~isempty(CLIM);
		CLIMlog = [-1 log10(CLIM(2))];
	else
		CLIMlog = CLIM;
	end
end

HDFfilename = specfilename;


%----------------  values for the GBS codes
% The +/- 5 window seems good for the new SG smoothing
% Used +/- 10 to evalauate contrast

jjj = 10; % X half-width of moving window 
iii = 10; % Y half-width of moving window

jjj = 5; % X half-width of moving window
iii = 5; % Y half-width of moving window

%jjj = 3; % X half-width of moving window
%iii = 3; % Y half-width of moving window

tbin = 1; % number of scans to bin together for 2-time calcs
tbin = 3; % number of scans to bin together for 2-time calcs

% used to 
% Express time axis in growth amount (ML) (will be used in user_input_guts_AREA (extra section)
delay = 3; % delay in frames of growth start after valve switch
%tamount = 4.5; % total ML of growth
tamount = 4.6; % total ML of growth  << this likely needs changing per scan

% read in the TwoTime ROI
%TTROIS_2017_04data   %<< presumably this would go back to being called in user space
%--------------------------------------------


% Now Call the meat
mocvdFLAG=1;

user_input_guts_AREA;




% write out the 

%DOCUMENTfigures = ...
%char(...
%	[' '],...
%	['Figure 1: [0  1 2] Single image Linear'],...
%	['Figure 2: [0 1 2] Single image Log int'],...
%	['Figure 3:  [3] Movie'],...
%	['Figure 4: [0 2 4] sum over X : Figure 5 sum over Y'],...
%	['    if particular point range chosen (POINTSUMS) then limited to those ranges'],...
%	['Figure 6: [0 2 4] sum over X and over all points : Figure 7 sum over Y and all points'],...
%	['Figure 8: [0 2 5] ROIs per scan variable (ROIS like counters)'],...
%	['Figure 9: [6] and above, Images (summed over points']);
%disp(DOCUMENTfigures);

%%% more documentation I don't want to lose
% Do not use x with this sset
%x [-4] Plot valves
%x [-3] Show raw images, harmonc, window and contamination, (summed over all points)
%x [-2] Show raw images, harmonc, window and contamination, integrated over YCOL 
%x [-1] Show raw images, harmonc, window and contamination, integrated over XROW
% [0] Just plot one SINGLE (spec point) image (Iwin) background subtracted and norm
% [1] Use as [0 1] to figure out ROIs to apply to image
% [2] Use as [0 2] to show the ROIS that are applied (ROIS below)
% [3] Make movie of Iwin (normed and log if chosen in flags sometimes you just
%			have to keep trying (close all figures) a few times to not get
%			the size error, it is inconsistent.
% [4] Use as [0 2 4] Image and line plots - Image summed over X and over Y  and 
%			also summed over the POINTSUMS if not empty)
% [5] Use as [0 2 5] ROIS (summed as counters, lines) VS Scan variable
% [6] Plots a 'single' image, but summed over %% POINTSUMS 
%		 It will make a new  figure for each POINTSUMS range
% [7] Plot the flatfield (show ROIS on it);
%x [7] Not implemented yet (2 time correlation)
%x [8] Look for static speckle after growth

% if any(WhatToDo==0)  % expect 1 file, plot the HDF image
%end
