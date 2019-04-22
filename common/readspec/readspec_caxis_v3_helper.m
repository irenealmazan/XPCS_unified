function [OM,scninfo,OMspecial] = readspec_caxis_v3_helper(scnheader,fileheader,comments);
% [OM,scninfo,OMspecial] = readspec_caxis_v3_helper(scnheader,fileheader,comments);
% version 20160213 
%
% OUTPUT (structures OM, scninfo, and OMspecial) 
%		OM  with following fields (shown with examples of what is in field) 
%		    angles_label: 'Nu  Chi  Eta  Delta'
%		    spectrometername: 'caxis_12id'
%	        cparam: [4.7580 4.7580 12.9900 90 90 120]
%		    h0: [3 0 0]
%		    h1: [1 1 3]
%		    a0: [18.4718 11.9530 -0.0988 0.0201]
%		    a1: [10.6740 39.4586 0.0900 5.8924]
%
%		OMspecial with following fields 
%           mode: 0
%			spectrometername: 'caxis_12id'
%           alphatarget: 0.12
%           betatarget: 1.2
%           sigma: -0.2404
%           tau: -37.4140
%           hklaz: [-0.0014 0.0011 1.0000]
%
%		scninfo   with the following fields
%			docu: [2x19 char]
%			spectrometername: 'caxis_12id'
%			mode: 0
%			filename: 'jul0915_2'
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
%			OM: << has the OM structure of above
%			Energy_i: 2.8000e+04
%			OMspecial: << has the OMspecial structure of above
%			datereduced (date this program was run to get this information)
%
% INPUT - uses last three outputs from readspecscan
%		[outdata, outpts, scantype, scandate, ncols, collabels, 
%			fileheader, scnheader, comments]   = readspecscan(name,scan)
%% sevchex specific comments 
		%% Should work for sevchex files (dec1314_1 and later in the Dec 2014 commisioning run) on 12id-d 
		%% Prior dec1314_1, sevchex files headers were a bit different

%% To create similar function specific for other spec files
	%%	readspecscan should work for all spec files
	%%	this helper function is designed so that only 2 subfunctions will need to be altered (or checked).


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

function [OM,GeoMot] = makeOM(G1)
	% *GeoMot* indexes where the geometry motors labeled by "angles_label" are located in the P0/P1/P2...scan header. See the #O label lines to determine the order of the motors. To specify indices, assume the P lines all concatenated as one long vector. 
	% *angledocu* - Make all functions consistent with using an order consistent with the printout of the  'pa' command. 
	% *angledocu* - REQUIRED, two spaces between each distinct motor name. (this is because too many systems use motors with two word names, e.g. "Two Theta". We cannot read minds, so we must assume that if there is only one space, the two words go together as one motor.
	
	% caxis Fileheader 2014  
	% #O0       Nu       Chi       Eta     Delta 
	% caxis Fileheader July 2015
	% #O0       Nu       Chi       Eta     Delta 

	GeoMot = [1:4]; 
	%%% MUST KEEP THE caxis motor name ORDER AS IN NAME LIST BELOW
	%%% IF THE name ORDER HAS CHANGED IN THE #00 line above
	%%% THEN RENUMBER the GeoMot vector to reflect new #O0 positions 
	%%        (e.g., GeoMot=[2 4 1 3] or what will find the motor names.
	%%
	%%% MUST KEEP THE caxis motor name ORDER AS IN NAME LIST BELOW
	OM.angles_label  = 'Nu  Chi  Eta  Delta' % 2 spaces between each separater motor name!
	OM.spectrometername = 'caxis_12id'; 
	% make these all row vectors for purposes of easy docu in structure
	% the cparams and h0, h1 positions are standard in the #G1 line
	OM.cparam 	= G1(1:6)';
	OM.h0		= G1(13:15)';
	OM.h1		= G1(16:18)';
	% a0 and a1 likely need changes to reflect number of geo motors.
	OM.a0		= G1([19:22])';
	OM.a1		= G1([25:28])';
	OM.a0E		= fhc./G1(31);  % energy that a0 was taken	
	OM.a1E		= fhc./G1(32);

end

%% OMspecial - make sure it has all the fields consistent with what is needed for the hkl2angle calculation functions for a particular spectrometer.

function [OMspecial] = extractspecial(scnheader,fileheader,scninfo)

	%index of lines
	G0ndx 	= intersect(strfind(scnheader(:,2)','G'), strfind(scnheader(:,3)','0')); 
	G4ndx 	= intersect(strfind(scnheader(:,2)','G'), strfind(scnheader(:,3)','4')); 
% #G4 line from caxis 2015
%#G4 3.000009054 -4.616810523e-06 0.04998622229 0.4428018736 0.1200424254
% 0.008755708022 30.524625 18.55257825 90.05638201 -0.39904 
% 38.12169 0.1200340486 14.42850009 0 0 0 0 -180 -180 -180 -180 0 1

% $G4 line from sevchex 2016
% #G4 -0.0002115385456 7.562283352e-05 0.5000106174 0.5165991452 1.509097462 
% 1.34938657 2.858502938 0.9999999949 0.15886 -51.3596 
% 0 0 0 0 0 0 -180 -180 -180 -180 -180 -180 -180
	% extract the numbers
	G0=sscanf(scnheader(G0ndx,[4:end]),'%g');
	G4=sscanf(scnheader(G4ndx,[4:end]),'%g');

	OMspecial.mode = round(G0(1));	
	OMspecial.alphatarget	= G4(5);
	OMspecial.betatarget	= G4(6);
	OMspecial.sigma		= G4(10);  
	OMspecial.tau		= G4(11); 
%	OMspecial.sigma		= G4(9);  
%	OMspecial.tau		= G4(10); 
	OMspecial.hklaz		= G0(4:6)';  % as row vector
% for purposes of using the sevc to calculate angles for a constant 'chi' scan, mode 4
%   (will need to use the 'chi' caxis value as 'phi' sevc
	OMspecial.Chi		= scninfo.geoangles_i(chan2col(scninfo.geoangles_label,'Chi'));
	% for purposes of using the sevc to calculate angles, convert Chi to Phi
	OMspecial.sign_nu = +1;	%%( hardwired, need to find which entry saves it	
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
%	G4	= sscanf(scnheader(G4ndx,[4:end]),'%g');  %%% save for later
		Pos = [];Olabels=[];
		for ii=[1:length(Pndx)];
			Posii = sscanf(scnheader(Pndx(ii),[4:end]),'%g');
			Olabelsii = fileheader(Ondx(ii),[4:end]);

			Pos = [Pos;Posii];
			Olabels = [Olabels,' ',Olabelsii];
		end
	
	QHKLi	= sscanf(scnheader(Qndx,[4:end]),'%g');

	[OM,GeoMot] = makeOM(G1);
	
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
	D = fix(clock)';
	plottime = ...
	[' [',datestr(now,1),'  ',datestr(now,15),'] '];
	DATEANDTIMESTRING = [plottime];
end

%%% example of output of caxis
%%% Unfortunately - caxis motor posiitons headers get clobbered randomly  and altered a lot by hasty user who don't think it important  to keep the GEOMETRY motors in the same place in the config.
%%% 1) Most important one is that the order and place of the important geometry motors changes in the config file when people move things without thinking of consequences
%%% 2) Also - from various macros interacting badly? and also the issues with missing command the last motor channels added in a config are NOT saved in the header until a freshstart is done.(Note all the 'Eta Eta Eta in the last #O7 row of this example). What this means is that there are some things we cannot retrieve from the file (depending on what motors or pseudo motors were stored in these regions of the config file.
%%
%%
% #F oct2714_1
% #E 1414413521
% #D Mon Oct 27 07:38:41 2014
% #C sput_141026a  User = bessrc
% #O0       Nu       Chi       Eta     Delta          Mu  TABLEVERT  TABLEHORZ
% #O1  DelS_TP   DelS_BM      mono  MONO_CHI    fpupvr    fpdnvr  MS0Bottom     MS0In
% #O2   MS0Out    MS0Top  MirBendDown  MirBendUp  MirVertUpOut  MirVertIn     Sam_Y     Sam_Z
% #O3    Sam_X  MirVertDownOut  MirTheta    MirChi   id12gap  Slit1Bot  Slit1Top  Slit1VGap
% #O4 Slit1VCen  Slit1out   Slit1in  Slit1HGap  Slit1HCen  Slit2Top  Slit2Bot  Slit2VGap
% #O5 Slit2VCen  Slit2InB  Slit2OutB  Slit2HGap  Slit2HCen    s0hgap    s0hcen    s0vgap
% #O6   s0vcen    MS1Bot    MS1Top   MS1VGap   MS1VCen    MS1InB   MS1OutB   MS1HGap
% #O7  MS1HCEN       Eta       Eta       Eta       Eta       Eta       Eta       Eta
% #C Book 92, Page 6
%
% #S 1  ascan  eta 0.829625 0.889625  30 1
% #D Mon Oct 27 07:41:26 2014
% #T 1  (Sec)
% #G0 2 0 1 0.003188262589 0.001725418314 0.9999934289 0 0 0 1
% #G1 3.905 3.905 3.905 90 90 90 1.609010322 1.609010322 1.609010322 90 90 90 2 0 1 2 0 3 22.851875 19.385 0.66775 11.2695 0 0 22.321875 43.479625 0.67375 35.881375 0 0 0.7749026676 0.7749026676
% #G3 1.600973911 -0.1606136007 -0.0001533800463 0.1606084002 1.600933065 -0.01150997611 0.001301551559 0.01143717793 1.608969146
% #G4 0.0001853693935 0.0004461859228 0.199991721 0.7749026676 1.152885409 1.121108988 -0.0006875 2.274000415 -96.87307652 0.3372 60.42435 0.9999968125 0 0 0 0 0 -180 -180 -180 -180 0 1
% #Q 0.000185369 0.000446186 0.199992
% #P0 0.001375 0 0.859625 2.274  0.00025 160.69996 73.399938
% #P1 3.299885 2.699935 7.0983204 0.1725 -3.934911 -5.065051 -0.6999998 -1.3999483
% #P2 1.400006 0.7 0.5008 0 0.50007 -1.8e-05 2.6000075 -2.111
% #P3 0.5000625 -0.50007 5.002 1.495 14.148786 -0.6825 1.7025 1.02
% #P4 -1.1925 4.18 5.8375 10.0175 0.8275 -0.0025 1.0025 1
% #P5 0.5025 1.8375 0.1625 2 -0.8375 2.0000004 2.5524764e-08 0.50000021
% #P6 0.028658151 9.475 1.625 11.3 0 444.51783 444.51783 444.51783
% #P7 444.51783 0.859625 0.859625 0.859625 0.859625 0.859625 0.859625 0.859625
% #CA Amp 1: 1.0e-06 (0.000e+00); Amp 2: 1.0e-06 (0.000e+00); Amp 3: 2.0e-08 (2.600e-10); 
% #N 62
% #L Eta  H  K  L  Epoch  Sec  ion1  Hubmon  scan_bar  Filters  Correct  Qmag  Energy  Graze  p_peak  p_total  p_del_s  p_nu_s  p_image  Vort1  Vort2  Vort3  Vort4  Vort5  Vort6  Vort7  Vort8  Vort9  Vort10  Vort11  Vort12  Vort13  Vort14  Vort15  Vort16  Vort17  Elastic  v_icr  v_ocr  v_real  v_live  v_dead  v_sum  p_xcen  p_xsig  p_ycen  p_ysig  power  tempsp  gun_on  setpw1  fwpw1  repw1  DCV1  setpw2  fwpw2  repw2  DCV2  Gamry_V  Gamry_I  ion0  p_point
% 0.829625 0.000174917 0.000341993 0.199992 167 1.0000006 471110 315525 0 2 504326.14 0.32178968 15.999971 1.1228855 327241 11910091 6687308 6400759 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 132135 0 0 0 0 55.1 750 553 0 0 0 0 0 0 0 0 0 0 392883 52615
% 0.831625 0.000175613 0.00034894 0.199992 168 1.0000006 469514 313326 0 2 527397.75 0.32178968 15.999971 1.1248855 343805 11818011 6632616 6349716 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 132135 0 0 0 0 55.1 750 553 0 0 0 0 0 0 0 0 0 0 391786 55022
%%%% etc

