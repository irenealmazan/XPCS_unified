function [OM,scninfo,OMspecial] = readspec_sevc_v1_helper(scnheader,fileheader,comments);
% designed to do various things to pull out information from headers and comments
% OUTPUT = readspec_sevc_v1_helper(filenames,scans,datapath,motorparams);
% Specific to the sevcfiles (dec1314_1 and later in the Dec 2014 commisioning run) on 12id-d 
% The few earlier files were a bit different as arranging motors, counters
% 2014-Dec-13 Carol Thompson cthompson@niu.edu

% Using output from the readspecscan
%[outdata,outpts,scantype,scandate,ncols,collabels,scnheader,fileheader,comments] = readspecscan(name,scan)

	% Extract sigma, tau, alpha, OM, EPOCH, motorpositions, motorlabels 

% OM, OMspecial
	[OM, scninfo] = extractOM(scnheader, fileheader);
	[OMspecial] = extractspecial(scnheader,fileheader);
	

end

%%%  HELPER FUNCTIONS %%%

function [OM,scninfo] = extractOM(scnheader,fileheader)
% The following will likely work for almost all more recent spec files
% and the only change will be to change the following line so that it specifies
% the crystallographic geometry motors and order that match the orientation matrix
% For example, in fourc, the motor names and order in the OM should be tth, th, chi, phi
% so if the motor slots show up as in line below, then GeoMot = [5 1 2 3]
%  #O0 Rho  Nu  Delta  Mu  Chi  Eta  Phi  DelRot
%  #O1 DelTrans
% If you are so unlucky that they are spread over #O0, #O1, #O2, then pretend all
% positions concatenated and use those numbers.

	GeoMot = [1:7];
	scninfo.docu  = 'Rho  Nu  Delta  Mu  Eta  Chi  Phi: sevchex'; 
	
%  G0 and G1 lines exist in scan headers in very old spec file (circa 2000) 
%	G0 is entries about the geometery, G1 is orientation matrix and energy when found
%	G3 (UB matrix, redundant information to G1) 
%	and G4 (more info about complex geo parameters and scan specific values
%	For per scan information in G3 and G4, we extract at end to build in some
%		sensible handling of old files

	%index of lines - note, older spec files will only have G0, and G1 
	G0ndx 	= intersect(strfind(scnheader(:,2)','G'),strfind(scnheader(:,3)','0'));
	G1ndx 	= intersect(strfind(scnheader(:,2)','G'),strfind(scnheader(:,3)','1')); 
%	G4ndx 	= intersect(strfind(scnheader(:,2)','G'),strfind(scnheader(:,3)','4'));
	Pndx	= strfind(scnheader(:,2)','P');
	Ondx	= strfind(fileheader(:,2)','O');
	Qndx	= strfind(scnheader(:,2)','Q');
	
	% extract the numbers
	G0	= sscanf(scnheader(G0ndx,[4:end]),'%g');
	G1	= sscanf(scnheader(G1ndx,[4:end]),'%g');
%	G4	= sscanf(scnheader(G4ndx,[4:end]),'%g');
		Pos = [];Olabels=[];
		for ii=[1:length(Pndx)];
			Posii = sscanf(scnheader(Pndx(ii),[4:end]),'%g');
			Olabelsii = fileheader(Ondx(ii),[4:end]);

			Pos = [Pos;Posii];
			Olabels = [Olabels,' ',Olabelsii];
		end
	QHKLi	= sscanf(scnheader(Qndx,[4:end]),'%g');

	% the orientation matrix information from the particular scan
	OM.cparam 	= G1(1:6);
	OM.h0		= G1(13:15);
	OM.h1		= G1(16:18);
	%angles, and energy at which they were found
	OM.a0		= {G1([19:24 33]),fhc./G1(31)};			
	OM.a1		= {G1([25:30 34]),fhc./G1(32)};
	OM.angledocu = '[Rho  Nu  Delta  Mu  Eta  Chi  Phi] : sevchex';

	scninfo.OM = OM;

% place holder for other modes, fixed beta, azimuth, theta, etc	% values going into scan (i.e., state before any action of scan) this is 'initial' values)
	scninfo.geomode	= round(G0(1));	   % should be integer
%	scninfo.Energy_i	= fhc./G4(4);	
	scninfo.HKL_i		= QHKLi(1:3);
	scninfo.geoangles_i	= [Pos(GeoMot)];
	scninfo.geoangles_label	= col2chan(Olabels,[GeoMot]);
	scninfo.allangles_i	= [Pos];
	scninfo.allangles_label	= col2chan(Olabels,[1:length(Pos)]);  %% 2016_02 column
		% Olabels has a lot of '' in it

% If older spec files - some parameters per scan were not kept (assumed more constant)
% and only G0 and G1 lines were defined
	G4ndx 	= intersect(strfind(scnheader(:,2)','G'),strfind(scnheader(:,3)','4'));

	if ~isempty(G4ndx)
		G4=sscanf(scnheader(G4ndx,[4:end]),'%g');
		scninfo.Energy_i	= fhc./G4(4);
	else
		scninfo.Energy_i = OM.a0{2}
		scninfo.docu = 'OLD spec file, assume energy at start of scan equivalent to energy during OM determination';
	end 	

	scninfo.docu = 'more info';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Separated these parameters which very specific to diffractometers or its modes
%	speciality of diffractometer or modes of this diffractometer. 
function [OMspecial] = extractspecial(scnheader,fileheader)

	%index of lines
	G0ndx 	= intersect(strfind(scnheader(:,2)','G'), strfind(scnheader(:,3)','0')); 
	G4ndx 	= intersect(strfind(scnheader(:,2)','G'), strfind(scnheader(:,3)','4')); 
	
	% extract the numbers
	G0=sscanf(scnheader(G0ndx,[4:end]),'%g');
	G4=sscanf(scnheader(G4ndx,[4:end]),'%g');

	OMspecial.alpha_i	= G4(5);
	OMspecial.beta_i	= G4(6);
	OMspecial.sigma		= G4(10);
	OMspecial.tau		= G4(11);
	OMspecial.hklaz		= G0(4:6);
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

function [colname,numofcols,channames_as_column] = col2chan(chan,colchoose)
% spec appears to put 2 spaces between each, with one space between #L and the start
% This will also accommodate ignorant people who named motors or detectors with spaces
% 'Two Theta'   (can be distinguished from
% 'Two  Theta' which is 'Two' and 'Theta'
% #L angle1  angle2  angle3 

	% 'Label1  Label2  Label Three' >> 'Label1','Label2','Label Three'
	chan_rephrased = [''''   (strrep(chan,'  ',''','''))    ''''];
		channames_as_columns = eval(['strvcat(',chan_rephrased,')']);
		numofcols = length(chan(:,1));   

	colname = strtrim(channames_as_columns(colchoose,:));
	

end

function [CHANCOLUMN,channelchoosematch] = chan2col(chan,channelchoose,quietflag)
% Columnnumber = chan2col(Speclabels,Chooselabel)
if nargin<3;quietflag=0;end

chan_rephrased = [''''   (strrep(chan,'  ', ''','''))    ''''];
	channames_as_column = strtrim(eval(['strvcat(',chan_rephrased,')']));     % this works OK for all labels
	numofcols = length(channames_as_column(:,1));   %just the total number of variables in the LABELS at start of scan	
	chan_rephrased2col = [ones(numofcols,1)*[''''], channames_as_column,  ones(numofcols,1)*[''';']];
	chan_rephrased2 = char(reshape(chan_rephrased2col',1,length(chan_rephrased2col(:))));
	channames_as_cellarray = strtrim(eval(['{',chan_rephrased2,'}']));

CHANCOLUMN = [];
channelchoosematch = [];jj=1;
for ii = 1:length(channelchoose(:,1))
	THE_MATCH 	= find( strcmp(strtrim(channelchoose(ii,:)),channames_as_cellarray));
		if ~isempty(THE_MATCH)
			channelchoosematch(jj) = ii;
			jj=jj+1;
		elseif quietflag ==0;  % if match is empty tell us unless quiet is on (not zero)
			disp([''''  strtrim(channelchoose(ii,:))  ''''  ' choice is not parsed within chan: ',chan]);
		end
	CHANCOLUMN = [CHANCOLUMN;THE_MATCH];
end
end
%%%%%%%%%%%%
% #F ./data/dec1214_1
% #E 1418411019
% #D Fri Dec 12 13:03:39 2014
% #C sevcmovcd  User = s12mocvd
%
% #O0 Phi  Chi  Rho  Eta  Nu  Delta  Mu  DelRot
% #O1 DelTrans  s6hgap  s6hcen  s6vgap  s6vcen  mono  monochi  mirvuo
% #O2 mirvdo  mirvi  mirbu  mirbd  ms0tp  ms0bm  ms0in  ms0out
% #o0 phi chi rho eta nu del mu delR
% #o1 delT s6hgap s6hcen s6vgap s6vcen mono monochi mirvuo
% #o2 mirvdo mirvi mirbu mirbd ms0tp ms0bm ms0in ms0out
%
% #J0 Seconds  ion0  ion1  pind  scanbar  filters  p_peak  p_total
% #J1 p_point  p_del_s  p_nu_s  p_image  correct  
% #j0 sec ion0 ion1 pind scanbar filters p_peak p_total
% #j1 p_point p_del_s p_nu_s p_image correct 
% 
%
% #S 1  ascan  phi 0.9 1.1  5 1
% #D Fri Dec 12 13:03:53 2014
% #T 1  (Seconds)
% #G0 0 0 0 0 0.01 1 750 0 1 0 0
% #G1 1.54 1.54 1.54 90 90 90 4.079990459 4.079990459 4.079990459 90 90 90 1 0 0 0 1 0 0 9.12 0 0 0 0 0 9.12 0 0 0 0 1.54 1.54 4.56 94.56
% #G3 4.079990459 -3.504879987e-15 -6.561207576e-16 2.848759229e-15 4.079990459 -2.498273628e-16 0 0 4.079990459
% #G4 0.4737451021 0.007078287183 0 1.033139392 -0.5642329723 0.5669537212 18.2894 0.5729386977 0 -90 0 0 0 0 0 0 -180 -180 -180 -180 -180 -180 -180
% #Q 0.473745 0.00707829 0
% #P0 1 0 0 0 18.2894 0 9.0007 2
% #P1 332.784 2.7773374 1.2682284 6.3574695 -1.2544829 9.4828829 -0.3 3.999928
% #P2 3.999928 4.000068 0 0.5008 1 -1.0000001 -1 1
% #N 18
% #L Phi  H  K  L  Epoch  Seconds  ion1  scanbar  filters  p_peak  p_total  p_point  p_del_s  p_nu_s  p_image  correct  ion0  pind
% #C Fri Dec 12 13:03:55 2014.  Scan aborted after 0 points.
%
% #S 2  ascan  phi 0.9 1.1  5
