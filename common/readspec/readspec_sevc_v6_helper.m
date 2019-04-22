function [OM,scninfo,OMspecial] = readspec_sevc_v6_helper(scnheader,fileheader,comments)
% [OM,scninfo,OMspecial] = readspec_caxis_v6_helper(scnheader,fileheader,comments);
% version 20160213 
% version 20160309 geomotor #P ordered automatically figured out
% version 20160331 OMspecial added detarm (arm length used in detector pseudomotor calc
%
%
% OUTPUT (structures OM, scninfo, and OMspecial) 
%		OM  with following fields (shown with examples of what is in field) 
%		    angles_label: 'Rho  Nu  Delta  Mu  Eta  Chi  Phi'
%		    spectrometername: 'sevchex'
%	        cparam: [3.1884 3.1884 5.1852 90 90 120]
%		    h0: [0 0 2]
%		    h1: [1 1 0]
%		    a0: [0.2 -0.1205 11.219 0 5.902 0.0281 5]
%		    a1: [0.2 18.606 0.1115 0 0.2265 0  -13.3143]
%			a0E: 28000	  {energy when a1 found}
%			a1E: 28000    {energy when a0 found}
%
%		OMspecial with following fields 
%           mode: 0
%			alphatarget: 0.12
%			betatarget:  1.2
%			sigma: 0.1379.82439
%			tau:  79.82439
%			hklaz: [0.001 0.0049 0.999987]
%			Rho: 0.2
%			Mu: 0
%			Phi: 7.4
%			Chi: 0.5
%			sign_nu: +1
%			detarm: 750			
%
%		scninfo   with the following fields
%			docu: [2x19 char]
%			spectrometername: 'sevchex'
%			mode: 0         { choices are 0 to 7 }
%			filename: '2016_03'
%			filedate: 'Thu Jul 09 18:22:20 2015'
%			epoch0: 1.4365e+09
%			scnnum: 1
%			scndate: 'Thu Jul 09 19:38:23 2015'
%			scntype: 'ascan  eta 3.3019 4.5019  60 0.5'
%			geoangles_i: [4x1 double]
%			geoangles_label: [4x5 char]
%			allangles_i: [67x1 double]
%			allangles_label: [67x14 char]
%			HKL_i: [3x1 double]
%			OM: << is the OM structure of above
%			Energy_i: 2.8000e+04
%			OMspecial: << is the OMspecial structure of above
%			datereduced (date this program was run to get this information)
%
% INPUT - uses last three outputs from readspecscan
%		[outdata, outpts, scantype, scandate, ncols, collabels, 
%			fileheader, scnheader, comments]   = readspecscan(name,scan)
%% sevchex specific comments 
		%% Should work for sevchex files (dec1314_1 and later in the Dec 2014 commisioning run) on 12id-d 
		%% Prior dec1314_1, sevchex files headers were a bit different ordering motors

%% To create similar function specific for other spec files
	%%	readspecscan should work for all spec files
	%%	this helper function is designed so that only 2 subfunctions will need to be altered (or checked).

% 2014-Dec-13 Carol Thompson cthompson@niu.edu
% 2016-Feb-12 C Thompson, making OMspecial more consisitent with what would be needed with the h2a_spectrometername (which will be different for each, so this one is specific for the sevchex
% 2016-Feb-13 C Thompson made a0 and a1 fields back to only the angles,
%	separate the energies that they were found at in a0E and a1E fields
%
% Put in current version (20160213) of chan2col function so that this function can stand alone.

	[OM, scninfo]	= extractOM(scnheader, fileheader);
	[OMspecial]		= extractspecial(scnheader,fileheader,scninfo);
	
	scninfo.OMspecial = OMspecial;
	scninfo.reduceddate = dateandtime;
	

end

%%%%%%%%%%%%%%%%%  HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% functions "makeOM"   and   "extractspecial" 
%% will likely need editing for new spectrometer versions.
%% Be aware of what OMspecial needs (to be passed to hkl2angle routines
%% for that particular spectrometer.
%% The rest of the functions should be OK without change even for 1990-mid2000 versions of spec. Not sure about 1980 or prior versions.

function [OM,GeoMot] = makeOM(G1,Olabels)
	% *GeoMot* indexes where the geometry motors are in "angles_label" 
	%	see the #P and the #O lines, and the #G1 line
	% *angledocu* - Make all functions consistent with using an order consistent with the printout of the  'pa' command. 
	% *angledocu* - REQUIRED, two spaces between each distinct motor name. (this is because too many systems use motors with two word names, e.g. "Two Theta". We cannot read minds, so we must assume that if there is only one space, the two words go together as one motor. We cannot read 3 word motors ('Two Baggins Theta')
	
	%%% We use the'pa' command to define the order of the geometry motor labels
	%%% We carry this through all related  programs (and users most likely 
	%%% would manually enter motor positions this way)
	%%% However - there is a potential issue that the 
	%%%% positions as saved in the OM ine (#G1) are ordered differently
	%%% This affects OM.a0 and OM.a1 (the angles) 
	%%%% If necessary, reorder these below, 

	OM.angles_label  = 'Rho  Nu  Delta  Mu  Eta  Chi  Phi'; % 2 spaces between each separater motor name!
	OM.spectrometername = 'sevchex'; 
	% make these all row vectors for purposes of easy docu in structure
	% the cparams and h0, h1 positions are standard in the #G1 line
	OM.cparam 	= G1(1:6)';
	OM.h0		= G1(13:15)';
	OM.h1		= G1(16:18)';
	% a0 and a1 likely need changes to reflect number of geo motors.
	% Whoever writes program should make sure order here matches labels
	OM.a0		= G1([19:24 33])';
	OM.a1		= G1([25:30 34])';
	OM.a0E		= fhc./G1(31);  % energy that a0 was taken	
	OM.a1E		= fhc./G1(32);
	
	% There is error that occurs occasionally in file header 
	% where some geo motors are repeated in place of motors in
	% last line. We do not want those doubled here. So, 
	% in the included chan2col function, we only output the '1st' instance
	GeoMot = chan2col(Olabels,OM.angles_label);

end

%% OMspecial - make sure it has all the fields consistent with what is needed for the hkl2angle calculation functions for a particular spectrometer.

function [OMspecial] = extractspecial(scnheader,fileheader,scninfo)

	%index of lines
	G0ndx 	= intersect(strfind(scnheader(:,2)','G'), strfind(scnheader(:,3)','0')); 
	G4ndx 	= intersect(strfind(scnheader(:,2)','G'), strfind(scnheader(:,3)','4')); 
	
	% extract the numbers
	G0=sscanf(scnheader(G0ndx,[4:end]),'%g');
	G4=sscanf(scnheader(G4ndx,[4:end]),'%g');

	OMspecial.mode = round(G0(1));	
	OMspecial.alphatarget	= G4(5);
	OMspecial.betatarget	= G4(6);
%	OMspecial.sigma		= G4(10); these seemed to change a bit in sevc until late 2015 
%	OMspecial.tau		= G4(11); 
	OMspecial.sigma		= G4(9);  
	OMspecial.tau		= G4(10); 
	OMspecial.hklaz		= G0(4:6)';  % as row vector  
	OMspecial.Rho		= scninfo.geoangles_i(chan2col(scninfo.geoangles_label,'Rho'));
	OMspecial.Mu		= scninfo.geoangles_i(chan2col(scninfo.geoangles_label,'Mu'));
	OMspecial.Phi		= scninfo.geoangles_i(chan2col(scninfo.geoangles_label,'Phi'));
	OMspecial.Chi		= scninfo.geoangles_i(chan2col(scninfo.geoangles_label,'Chi'));
	OMspecial.sign_nu 	= G0(9);
	OMspecial.detarm	= G0(7);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The following should work without alteration %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [OM,scninfo] = extractOM(scnheader,fileheader)
%	Only #G0 and #G1 lines exist in very old spec file (circa 2000) 
%	G0 documents details about geometery, G1 is orientation matrix 
%	G3 (UB matrix, redundant information to G1) 
%	G4 (more info about geometry specific parameters and scan specific values
%	For information in G3 and G4, we use the information when it exists,
%	but have built in some handling for old files. 
% 	scninfo has been written to be compatible with all spec files even 1998

	%index of lines - note, older spec files will only have G0, and G1 
	G0ndx 	= intersect(strfind(scnheader(:,2)','G'),strfind(scnheader(:,3)','0'));
	G1ndx 	= intersect(strfind(scnheader(:,2)','G'),strfind(scnheader(:,3)','1')); 
%	G4ndx 	= intersect(strfind(scnheader(:,2)','G'),strfind(scnheader(:,3)','4'));
	Pndx	= strfind(scnheader(:,2)','P');
	Ondx	= strfind(fileheader(:,2)','O');
	Qndx	= strfind(scnheader(:,2)','Q');
	Sndx	= strfind(scnheader(:,2)','S');
		STndx = strfind(scnheader(Sndx,:),' ');
	SDndx	= strfind(scnheader(:,2)','D');
	FDndx	= strfind(fileheader(:,2)','D');
	Fndx	= strfind(fileheader(:,2)','F');
	Epndx	= strfind(fileheader(:,2)','E');
	
	% extract the numbers
	G0	= sscanf(scnheader(G0ndx,[4:end]),'%g');
	G1	= sscanf(scnheader(G1ndx,[4:end]),'%g');
%	G4	= sscanf(scnheader(G4ndx,[4:end]),'%g');  %%% 
		Pos = [];Olabels=[];
		for ii=[1:length(Pndx)];
			Posii = sscanf(scnheader(Pndx(ii),[4:end]),'%g');
			Olabelsii = fileheader(Ondx(ii),[4:end]);

			Pos = [Pos;Posii];
			Olabels = [Olabels,' ',Olabelsii];
		end
	
	QHKLi	= sscanf(scnheader(Qndx,[4:end]),'%g');

	[OM,GeoMot] = makeOM(G1,Olabels);
	scninfo.docu  = OM.angles_label;
	scninfo.spectrometername = OM.spectrometername;	



% 'initial' values at start of scan or file
	scninfo.mode	= round(G0(1));	   %% should be treated as integer
%	scninfo.Energy_i	= fhc./G4(4);  %% save for later, no #G4 old files
	scninfo.filename = strtrim(fileheader(Fndx,[4:end]));  %% within file, 
	scninfo.filedate = strtrim(fileheader(FDndx,[4:end]));
	scninfo.epoch0 = sscanf(fileheader(Epndx,[4:end]),'%g');
	scninfo.scnnum = sscanf(scnheader(Sndx,[4:end]),'%g');
	scninfo.scndate = strtrim(scnheader(SDndx,[4:end]));
	scninfo.scntype = strtrim(scnheader(Sndx,[STndx(3)+1:end]));
%	scninfo.scntype = scntype;
	scninfo.geoangles_i	= [Pos(GeoMot)];
	scninfo.geoangles_label	= col2chan(Olabels,[GeoMot]);
	scninfo.allangles_i	= [Pos];
	scninfo.allangles_label	= col2chan(Olabels,[1:length(Pos)]);  
	scninfo.HKL_i		= QHKLi(1:3);
	scninfo.OM = OM;

% Information in #G4 only exists in later spec files.
	G4ndx 	= intersect(strfind(scnheader(:,2)','G'),strfind(scnheader(:,3)','4'));

	if ~isempty(G4ndx)
		G4=sscanf(scnheader(G4ndx,[4:end]),'%g');
		scninfo.Energy_i	= fhc./G4(4);
	else
		scninfo.Energy_i = OM.a0E
		scninfo.docu = char(scninfo.docu,'OLD spec file, assume energy at start of scan is equivalent to energy during OM determination');
	end 	

	scninfo.docu = char(scninfo.docu,'more info');
end


function hc = fhc
% output  12398.419292;  %hc in eV-angstroms best from CODATA constants;
%  01-March-2013 Carol Thompson   cthompson@niu.edu
hc = 12398.419292;
	% CODATA http://physics.nist.gov/cuu/Constants/index.html
	% CODATA Internationally recommended values of the constants (2010)
	% h = 6.62606957e-34;	%Joule-sec
	% c  = 299792458;	% m/sec
	% e = 1.602176565e-19;	% Coulomb/e
	% hc  = h*c/e;		% eV.m  (J/e(Coulombs) = eV)
	% hc  = h*c*1e10/e	% eV.Angstrom  (1e10 Ang / meter)
	% hc = 12398.4192920042	% eV.Angstrom
end

function [CHANCOLUMN,channelchoosematch] = chan2col(chan,channelchoose,quietflag)
% Columnnumber = chan2col(Speclabels,Chooselabel)
% from CT external version 20160213  (to make this function self-contained)
% see the chan2col function for more documentation
% Note extra line added that is not in version 20160213 so that this version of 
% chan2col will only give the first instance of a match, and not all 
% (needed to ignore header issues)

if nargin<3;quietflag=0;end

% replace double spaces with commas and  turn into column of strings
if length(chan(:,1))==1; 
	chan_rephrased = [''''   (strrep(chan,'  ', ''','''))    ''''];
	channames_as_column = strtrim(eval(['strvcat(',chan_rephrased,')']));
else  % or it already is column array
	channames_as_column = chan;
end
    % Then turn the column into a cell array
	numofcols = length(channames_as_column(:,1));   
	chan_rephrased2col = [ones(numofcols,1)*[''''], channames_as_column,  ones(numofcols,1)*[''';']];
	chan_rephrased2 = char(reshape(chan_rephrased2col',1,length(chan_rephrased2col(:))));
	channames_as_cellarray = strtrim(eval(['{',chan_rephrased2,'}']));
	
if length(channelchoose(:,1))==1;
	chanchoose_rephrased = [''''   (strrep(channelchoose,'  ', ''','''))    ''''];
	chanchoosenames_as_column = eval(['strvcat(',chanchoose_rephrased,')']);
	channelchoose = chanchoosenames_as_column;
end

CHANCOLUMN = [];
channelchoosematch = [];jj=1;
for ii = 1:length(channelchoose(:,1))
	THE_MATCH 	= find( strcmp(strtrim(channelchoose(ii,:)),channames_as_cellarray));
		if ~isempty(THE_MATCH)
			THE_MATCH = THE_MATCH(1);  %%% CT added specifically to kick out multiples
			channelchoosematch(jj) = ii;
			jj=jj+1;
		elseif quietflag ==0;  % if match is empty tell us
			disp([''''  strtrim(channelchoose(ii,:))  ''''  ' choice is not parsed within chan: ',chan]);
		end
	CHANCOLUMN = [CHANCOLUMN;THE_MATCH];
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DATEANDTIMESTRING = dateandtime
%function DATEANDTIMESTRING = dateandtime
% function that outputs date and time in string
	D = (clock);
	DATEANDTIMESTRING = ...
	['[',datestr(D,1),'  ',datestr(D,15),']'];
end

%% Example of file and scan header used for this sevc.

% #F 2016_0208_3
% #E 1454982072
% #D Mon Feb 08 19:41:12 2016
% #C pddewet1  User = mocvd
% 
% #O0 Rho  Nu  Delta  Mu  Eta  Chi  Phi  DelRot
% #O1 DelTrans  s6hgap  s6hcen  s6vgap  s6vcen  mono  monochi  mirvuo
% #O2 mirvdo  mirvi  mirbu  mirbd  ms0tp  ms0bm  ms0in  ms0out
% #O3 inczui  inczuo  inczd  TX  TY  TZ  dettran  frontzu
% #O4 frontzd  chambz1  chambz2  chambz3  mostab  
% #o0 rho nu del mu eta chi phi delR
% #o1 delT s6hgap s6hcen s6vgap s6vcen mono monochi mirvuo
% #o2 mirvdo mirvi mirbu mirbd ms0tp ms0bm ms0in ms0out
% #o3 inczui inczuo inczd TX TY TZ dettran frontzu
% #o4 frontzd chambz1 chambz2 chambz3 mostab 

% #J0 Seconds  ion0  ion1  ion2  hubmon  hexmon  laser  tempC
% #J1 filters  correct  p_point  p_image  p_total  p_peak  p_small  p_med
% #J2 p_lar  p_xcen  p_ycen  p_xsig  p_ysig  valve_p  valve_v  mirbd_s
% #J3 mirbu_s  
% #j0 sec ion0 ion1 ion2 hubmon hexmon laser tempC
% #j1 filters correct p_point p_image p_total p_peak p_small p_med
% #j2 p_lar p_xcen p_ycen p_xsig p_ysig valve_p valve_v mirbd_s
% #j3 mirbu_s 


% #S 1  ascan  TZ -1.3766 -0.9766  50 0.5
% #D Mon Feb 08 19:50:19 2016
% #T 0.5  (Seconds)
% #G0 3 0 0 0.001136619605 0.004916311814 0.9999872689 750 0 1 1 1
% #G1 3.1884 3.1884 5.1852 90 90 120 2.275497885 2.275497885 1.211753704 90 90 60 0 0 2 1 1 0 0.2 -0.1205 11.219 0 5.902 0.0281 0.2 18.606 0.1115 0 0.2265 0 0.5165996245 0.5165996245 5 -13.3143
% #G3 1.381509538 2.256476873 -0.01382933501 -1.807934276 0.2926173947 0.00710910454 0.02637564337 0.02403766141 1.211653932
% #G4 0 0 0 0.5165996245 -0.1069842966 0.1069842966 0 1 -0.31189 79.82439 0 0 0 0 0 0 -180 -180 -180 -180 -180 -180 -180
% #Q 0 0 0
% #P0 0.2 0 -0.2 0 0 0 0 -0.2
% #P1 -2.618 0.29998746 -3.8617886e-08 0.029979632 -5.2845529e-09 4.7254531 -1.035 0.548418
% #P2 -0.345862 0.099224 -2.1904 -12 0.5 -0.50000002 -0.5 0.5
% #P3 -13.1699 -13.1699 -14.37005 0 0.2601 -1.1766 -75 -3.729991
% #P4 -4.986655 -16.97 -16.97 -16.97 23409 

% #G0 line for sevc (as of 2015-05) - starting with 0 as index (not 1)
% /* Geometry Parameters */
% #define f_mode          (*gparam[0])    /* mode    spectrometer mode */
% #define f_sector        (*gparam[1])    /* sector  (not used) */
% #define f_freeze        (*gparam[2])    /* frz     frozen flag */
%
% /* Azimuthal reference (surface normal in hkl) */
% #define f_azref0        (*gparam[3])    /* h_az */
% #define f_azref1        (*gparam[4])    /* k_az */
% #define f_azref2        (*gparam[5])    /* l_az */
%
% #define f_len           (*gparam[6])    /* g_len  del arm length */
% #define f_zrot_p        (*gparam[7])    /* g_zrot true if CCW z rotations positive */
% /* Note: internal equations here are based on CW z rotations (Nu, Mu, Phi) */
% /* Signs are flipped on input and output from c_angle for CCW (if f_zrot_p true) */
% /* f_zrot_p should be false for MOCVD spectrometer, true for front spectrometer */
% #define f_NU_p          (*gparam[8])    /* g_NU true to choose positive c_angle[Nu] values */
% #define f_track         (*gparam[9])    /* g_track true when delR, delT tracked with del */
% #define f_sigtau        (*gparam[10])   /* set if sigtau set before azref */
