function [OM,scninfo,OMspecial] = readspec_sevc_v2_helper(scnheader,fileheader,comments)
% Pull out information from headers and comments (as read out using readspecscan
% [OM,scninfo,OMspecial] = readspec_sevc_v2_helper(scnheader,fileheader,comments);
%
% Expected to be using output from the readspecscan
% [outdata,outpts,scantype,scandate,ncols,collabels,scnheader,fileheader,comments]...
%			= readspecscan(name,scan)
%
% Specific to the sevcfiles (dec1314_1 and later in the Dec 2014 commisioning run) on 12id-d 
% The few earlier files were a bit different as arranging motors, counters
% 2014-Dec-13 Carol Thompson cthompson@niu.edu
%
% 2015-Dec-30 Carol Thompson  fixed sevc geo positions
%
% OUTPUT (structures)   (read from scan header or file header
%	OM, scninfo should have similar fields for all spectrometers 
%	OMspecial will have information specific to spectrometer
%
%	OM.cparam     (vector [a b c alpha beta gamma] (with angles in degrees)
%	OM.h0  and OM.h1  (are vectors [h k l] of each of the orientation (or0, or1)		
%	OM.a0  and OM.a1  (are 2 element cells {angles, Energy} 
%						where angles is vector of orientation angle 
%						(put in order as output by 'pa' command)
%						Energy is the energy when or0 (or1) input
%	OM.anglesdoeu		'TwoTheta  Theta  :  twocircle';  % angles, 
%								2 spaces between, and : spectrometer
%  
%		% _i refers to values prior to start of scan (in scan header)
%	scninfo.geomode		(0, 1, etc, most geometries have more than one mode
%	scninfo.HKL_i 		(HKL prior to scan
%	scninfo.geoangles_i (Positions of geometry relevant motors prior to scan
% 	scninfo.geoangles_labels  (Names of the geometry motors)
%	scninfo.allangles_i (Positions of all motors prior to scan)
%	scninfo.allangles_labels (Names of all motors - 2 spaces between each)
%	scninfo.Energy_i    (Energy prior to scan)
%	scninfo.OM   		(scaninfo.OM = OM from above)
%	scninfo.docu		(information)
%			
%   These are specific to sevchex although similar for most z-axis related
%	OMspecial.alpha_i	alphatarget prior to scan
%	OMspecial.beta_i	betatarget prior to scan
%	OMspecial.sigma		sigma 
%	OMspecial.tau		tau 
%	OMspecial.hklaz		sigmatau expressed as hkl azimuth (hkl normal to phi circle)
%	OMspecial.docu		usually spectometer, geometry

% Expected to be using output from the readspecscan
% [outdata,outpts,scantype,scandate,ncols,collabels,scnheader,fileheader,comments]...
%			= readspecscan(name,scan)

% Easy to find place for version info

	DOCU = 'sevc_v3 2016_02 CT'; 	

% Extract sigma, tau, alpha, OM, EPOCH, motorpositions, motorlabels 	
% OM, OMspecial
	[OM, scninfo] = extractOM(scnheader, fileheader);
	[OMspecial] = extractspecial(scnheader,fileheader);

	scninfo.docu = strvcat(DOCU,scninfo.docu);

end

%%%  HELPER FUNCTIONS %%%

function [OM,scninfo] = extractOM(scnheader,fileheader)
% The following will likely work for almost all more recent spec files
% and the only change will be to change the three line so that it specifies
% the crystallographic geometry motors and order that match the orientation matrix
% These iines are marked with %%%%% CHANGE %%%%%%%%%%
% For example, in fourc, the motor names and order in the OM should be TwoTheta, Theta, Chi, PHI
% so if the motor slots show up as in line below, then GeoMot = [5 1 2 3]
%  #O0 Theta  Chi  Phi  Sample  Theta  Tablex Tabley Tablez
%  #O1 Monochi  Mono ...
%
% If you are so unlucky that they are spread over #O0, #O1, #O2, then pretend all
% positions concatenated and use those numbers.
% Note - there are two spaces in between each motor name (to not be confused with single motor name with space

%%%%% CHANGE %%%%%%%%%%
	GeoMot = [1:7];
	scninfo.docu  = '[Rho  Nu  Delta  Mu  Eta  Chi  Phi] : sevchex'; 
%%%%%%%%%%%%%%%%%%%%%%

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
	%%%%% CHANGE %%%%%%%%%%
	%OM.a0		= {G1(19:22),fhc./G1(31)};			
	%OM.a1		= {G1(25:28),fhc./G1(32)};
	OM.a0		= {G1([19:24 33]),fhc./G1(31)};			
	OM.a1		= {G1([25:30 34]),fhc./G1(32)};
	OM.angledocu = '[Rho  Nu  Delta  Mu  Eta  Chi  Phi] : sevchex';
	%%%%%%%%%%%%%%%%%%%%%%%%%

	scninfo.OM = OM;

% place holder for other modes, fixed beta, azimuth, theta, etc	% values going into scan (i.e., state before any action of scan) this is 'initial' values)
	scninfo.geomode	= round(G0(1));	   % should be integer
%	scninfo.Energy_i	= fhc./G4(4);	
	scninfo.HKL_i		= QHKLi(1:3);
	scninfo.geoangles_i	= [Pos(GeoMot)];
	scninfo.geoangles_label	= col2chan(Olabels,[GeoMot]);
	scninfo.allangles_i	= [Pos];
	scninfo.allangles_label	= Olabels;

% If older spec files - some parameters per scan were not kept (assumed more constant)
% and only G0 and G1 lines were defined
	G4ndx 	= intersect(strfind(scnheader(:,2)','G'),strfind(scnheader(:,3)','4'));

	if ~isempty(G4ndx)
		G4=sscanf(scnheader(G4ndx,[4:end]),'%g');
		scninfo.Energy_i	= fhc./G4(4);
		scninfo.docu = strvcat(scninfo.docu,...
			'Energy of scan may or may not be equivalent to energy during OM determination');
	else
		scninfo.Energy_i = OM.a0{2}
		scninfo.docu = strvcat(scninfo.docu,...
			'OLD spec file, assume energy at start of scan equivalent to energy during OM determination');
	end 	

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

	%% note that scninfo.geomode will document the mode for the scan, 
	%% at minimum, these will be the minimium (per scan) to calculate angles
	
	%%%%% CHANGE AS NEEDED HERE FOR WHAT IS SPECIFIC DIFFFERENT SPECTROMETER %%%%%%%%%%
	OMspecial.alpha_i	= G4(5);
	OMspecial.beta_i	= G4(6);
	OMspecial.sigma		= G4(10);
	OMspecial.tau		= G4(11);
	OMspecial.hklaz		= G0(4:6);
	OMspecial.docu		= 'sevchex  12id-d APS';
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% #G1 1.54 1.54 1.54 90 90 90 4.079990459 4.079990459 4.079990459 90 90 90 1 0 0 0 1 0 0 9.12 0 0 0 0 0 9.12 0 0 0 0 1.54 1.54 4.56 94.56	%angles, and energy at which they were found

% #G1 
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
%
%
% #S 5  ascan  phi -14.021 -13.971  40 1
% #D Thu Dec 17 12:12:56 2015
% #T 1  (Seconds)
% #G0 4 0 0 -0.0006305432812 0.004509907608 0.9999896315 750 0 1 1 1
% #G1 3.186 3.186 5.178 90 90 120 2.277212008 2.277212008 1.213438646 90 90 60 0 0 2 1 1 0 0.2 -0.1227 11.212 4.999999987e-06 5.6071 0 0.2 18.6406 0 0 0 0 0.5165991452 0.5165991452 6.7 -13.9807
% #G3 1.36087797 2.261525422 -0.01297878623 -1.825778157 0.2657416868 0.0006173342611 0.01548552706 0.02405517342 1.213369077
% #G4 0.9992057545 0.9987577263 0.03825184812 0.5165991452 0.1859846882 0.1496829069 18.64516387 0.8978342175 -0.22164 39.20359 0 0 0 0 0 0 -180 -180 -180 -180 -180 -180 -180
% #Q 0.999206 0.998758 0.0382518
% #P0 0.2 18.6409 0.2005 0 0.0804 0 -13.991 0.2005
% #P1 2.6245 0.0029848171 -3.0487805e-08 0.0029852693 -3.4146341e-08 4.7254531 -0.9125 0.444928
% #P2 -0.444928 0 -13.9984 -3.8736 1 -1 -2 2
% #P3 -13.71981 -13.71981 -14.91996 0 -0.25 0.934 -7.4255 -3.729991
% #P4 -4.986655 -16.97 -16.97 -16.97 
% #U0 6553.000000 9175.000000 0.000000 3604.000000 0.000000 1638.000000 8847.000000 0.000000
% #U1 8847.000000 7249.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
% #U2 0.000000 0.000000 0.000000 0.000000 3276.000000 655.000000 2730.000000 2730.000000
% #U3 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
% #U4 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
% #U5 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
% #U6 0.000000 0.000000 31129.000000 0.000000 0.000000 0.000000 32768.000000 32768.000000
% #U7 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
% #U8 13094.000000 18355.000000 146.000000 7207.000000 0.000000 3266.000000 17702.000000 830.000000
% #U9 17688.000000 14488.000000 77.000000 7837.000000 0.000000 0.000000 0.000000 0.000000
% #U10 0.000000 0.000000 0.000000 0.000000 6544.000000 1301.000000 5451.000000 5455.000000
% #U11 12.000000 4.000000 0.000000 0.000000 266.000000 216.000000 201.000000 165.000000
% #U12 259.000000 279.000000 276.000000 10309.000000 1750.000000 31129.000000 1457.000000 8231.000000
% #U13 2131.000000 33167.000000 0.000000 5819.000000 5817.000000 0.000000 13115.000000 0.000000
% #U14 247.000000 253.000000 107.000000 264.000000 267.000000 269.000000 4990.000000 450.000000
% #U15 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
% #U16 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 3928.000000 1638.000000
% #U17 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
% #U18 0.000000 0.000000 6157.000000 0.000000 0.000000 0.000000 0.000000 0.000000
% #U19 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
% #U20 4095 1664 16896 64 64505 1 3840 8192
% #U21 61440 63871 33023 8127 1028 8192 341 16382
% #U22 1 68
% #U23 0 0 0 2
% #N 28
% #L Phi  H  K  L  Epoch  Seconds  ion0  ion1  hubmon  ion2  laser  tempC  filters  correct  p_point  p_image  p_total  p_peak  p_small  p_med  p_xcen  p_ycen  p_xsig  p_ysig  valve_p  valve_v  hexmon  p_lar

