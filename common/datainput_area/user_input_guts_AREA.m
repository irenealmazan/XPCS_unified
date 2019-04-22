% called by datainput_AREA which set ups various initializations, and also uses a sample file
% took out fixlogimage since matlab past 2016? will not croak in an image or pcolor
% also took out fixing NaN since now the sum and mean allow a 'omitnan' flag (that outputs similar to my outputs)
%		Note  - may need to check - not sure if I have currently got issues with NaN or zer
% 	if one of the pixels is NaN (or 0) (and it seems to do it sensibly by ignoring them)

% If ImageJ = 1; (then we are using ImageJ and spec indexing in ROI's and scan points (start at 0)
%			= 0; Then we are using matlab indexing in ROI's and scan points
ImageJ = 1; % true - using imageJ/spec indexing [start at 0] 


if ~ReUseFlag  % read in from scratch the sdata and idata

	function_init_readinnew;

else  % use the prior sdata and idata structures

	function_init_readinold;

end

	function_init_common



%%%%  CT took out the minus flags choices(used for the pixirad and windowing) %%%
%%%% and put in them in a separate script, I will comment it out as we are not using it
%%%% for the medipix anyway but it helps clean up this frankenstein file

%user_input_guts_MPX3_minusflags

%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(WhatToDo==0); % plot a 'single' image from the scan, typically used to check ROI's and such

	function0
	
end	


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(WhatToDo==1) %& ~isempty(ROIS) %% open up make ROIS
	myROIS = makerois(4,gcf);
	xlabel(XCOL);ylabel(YROW);  % must have done the WhattoDo =0
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(WhatToDo==2)  % plot the ROIS on figure (above, (from 1) and also make a cleam figure showing them

	function2
		
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(WhatToDo==3)  %% make movie  (frame must be 560x413?)

	function3
	
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(WhatToDo==4) %% make summed images 

	function4  % this is very long one, still script
	
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% note if the ROI has NaN in any of the pixels, this will mess up and give NaN
%%% need to fix sumrois to fix that
if any(WhatToDo==5)  %% make summed lines of ROI's on images

	function5
	
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(WhatToDo==6);   %to make it automatic for summs   % ? is Xsteps missing one point

	function6;  % in prep for making function (currently script)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(WhatToDo==7);  % plot flatfield FLATFIELDFLAG

	% currently just script - getting ready to turn into functions
	function7
 
end	

%---------------------  GBS sections WhattoDo 8 and above
 % user_input_guts_2timesection
%---------------------------
% set up the idata to pass out

% Create Unique ID so that we can associate that the sdata and the idata have not been mixed
	rng('shuffle');ID = rand(5,5);
    Unique_ID = (uint64(floor(rand.*1e18)));   
    OurDate	  = now;

sdata.Uid		= Unique_ID;   % should check this when reusing sdata to work with altering idata 
sdata.ReduceDate = datestr(OurDate);

idata.IInb 		= IInormb;
idata.XCol 		= XCOLpts;
idata.YRow 		= YROWpts;
idata.Zt		= timestampX;
idata.Z			= SPECpts;
idata.XColstr	= XCOL;
idata.YRowstr	= YROW;
idata.Normalization 	= Norm;
idata.Normalizationstr	= NORMdoc;
idata.Uid		= Unique_ID;
idata.ReduceDate = datestr(OurDate);
idata.FLATFIELDFLAG = [FLATFIELDFLAG];
idata.FLATFIELD = FLATFIELD;






