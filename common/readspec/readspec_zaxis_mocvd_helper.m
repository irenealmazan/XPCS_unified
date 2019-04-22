function OUTPUT = readspec_zaxis_mocvd_helper(scnheader,fileheader,comments);
% designed to do various things to pull out information from headers and comments
% OUTPUT = readspec_zaxis_mocvd_helper(filenames,scans,datapath,motorparams);
% Specific to the zaxis files (post 2002 and prior to final run July 2013) on 12id-d 
% However, it will likely work for the  huber at 12id set up running 'caxis' program (since 200X)
%	as the headers are set up the same and only four angles usedge (Nu Chi Eta Delta)
%		([Nu and Delta] is similar to [Delta and Gamma] in zaxis, and [Chi and Eta] to  [Theta and Mu]
% 2013-Jul-30 Carol Thompson cthompson@niu.edu

% Using output from the readspecscan
%[outdata,outpts,scantype,scandate,ncols,collabels,scnheader,fileheader,comments] = readspecscan(name,scan)

	% Extract sigma, tau, alpha, OM, EPOCH, motorpositions, motorlabels 

% OM, OMspecial
	[OM, scninfo] = extractOM(scnheader, fileheader);
	[OMspecial, scninfo] = extractspecial(scnheader,fileheader);
	

end

%%%  HELPER FUNCTIONS %%%

function [OM,scninfo] = extractOM(scnheader,fileheader)
% The following will likely work for almost all more recent spec files
% and the only change will be to change the following line so that it specifies
% the crystallographic geometry motors and order that match the orientation matrix
% For example, in fourc, the motor names and order in the OM should be tth, th, chi, phi
% so if the motor slots show up as in line below, then GeoMot = [5 1 2 3]
% #O0    th       chi       phi       t3x       tth    analth   analtth
% If you are so unlucky that they are spread over #O0, #O1, #O2, then pretend all
% positions concatenated and use those numbers.

	GeoMot = [1:4];

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
			Posii = sscanf(scnheader(Pndx(ii),[4:end]),'%g')
			Olabelsii = fileheader(Ondx(ii),[4:end])

			Pos = [Pos;Posii];
			Olabels = [Olabels,Olabelsii];
		end
	QHKLi	= sscanf(scnheader(Qndx,[4:end]),'%g');

	% the orientation matrix information from the particular scan
	OM.cparam 	= G1(1:6);
	OM.h0		= G1(13:15);
	OM.h1		= G1(16:18);
	%angles, and energy at which they were found
	OM.a0		= {G1(19:22),fhc./G1(31)};			
	OM.a1		= {G1(25:28),fhc./G1(32)};

	scninfo.OM = OM;

% place holder for other modes, fixed beta, azimuth, theta, etc	% values going into scan (i.e., state before any action of scan) this is 'initial' values)
	scninfo.geomode	= round(G0(1));	   % should be integer
%	scninfo.Energy_i	= fhc./G4(4);	
	scninfo.HKL_i		= QHKLi(1:3);
	scninfo.geoangles_i	= [Pos(GeoMot)];
	scninfo.geoangles_label	= col2chan(OLabels,[GeoMot]);
	scninfo.allangles_i	= [Pos];
	scninfo.allangles_label	= OLabels;

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

	scninfo.docu = 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Separated these parameters which very specific to diffractometers or its modes
%	speciality of diffractometer or modes of this diffractometer.
function [OMspecial] = extractspecial(scnheader)

	%index of lines
	G0ndx 	= intersect(strfind(scnheader(:,2)','G', strfind(scnheader(:,3)','0')); 
	G4ndx 	= intersect(strfind(scnheader(:,2)','G', strfind(scnheader(:,3)','4')); 
	
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
%%%%%%%%%%%%
% Example of header structures for the zaxis_mocvd station (post 2002)
% #F jul0813_2
% #E 1373327025
% #D Mon Jul 08 18:43:45 2013
% #C n130704a  User = mocvd
% #O0    Delta     Theta        Mu     Gamma   GamTran    GamRot        sh        sv
% #O1    tabhu     tabhd    tabvui    tabvuo     tabvd     mirbd     mirbu    mirvuo
% #O2    mirvi    mirvdo      mono   monochi     s0bot     s0top     s0out      s0in
% #O3     s8vt      s8vb     s8vap     s8vce      s5vt      s5vb     s5vap     s5vce
% #O4    s6hce     s6hap     s6vce     s6vap     flhap     flhce      flvb      flvt
% #O5    flvap     flvce      s1vb      s1vt     s1vap     s1vce      s1hl      s1hr
% #O6    s1hap     piezo    fpdnvr    fpupvr     s1hce  
% #C Book 184, Page 119
%
% #S 1  ascan  th 9.543 10.343  20 1
% #D Mon Jul 08 18:44:21 2013
% #T 1  (Seconds)
% #G0 0 0 0 0.9998358466 0.005865812378 -0.01714269856 426.17 0 1 1
% #G1 3.186 3.186 5.178 90 90 120 2.277212008 2.277212008 1.213438646 90 90 60 1 0 -4 2 0 -4 19.4142 16.8378 0.27844 8.993214286 0 0 19.573 23.4474 0.28159 18.41142857 0 0 0.438105327 0.438105327
% #G3 -0.02073676846 -0.1863660923 -1.208546398 -0.006670972784 1.960914616 -0.1082623778 2.277107818 1.142705594 -0.01132294308
% #G4 0.5042600719 -0.0193806687 -1.987199137 0.438105327 0.06178727928 4.529036538 5.1349 10.64878069 87.51540034 0.1701 -99.41293 0.2378390714 0 0 0 0 0 -180 -180 -180 -180 0
% #Q 0.50426 -0.0193807 -1.9872
% #P0 9.6162 9.943 0.23188 4.3610714 32.499375 4.361 4.4196088e-08 0
% #P1 3.574603 3.574603 -3.014599 -3.0146078 -3.0146062 6.3536 -4.4496 0.443032
% #P2 0 -0.443032 4.006175 0.7525 0.79999999 1.2 1.2 0.9
% #P3 0.3825 0.3625 0.49893 -0.8675 1.2625 3.0975 4.36 0.9175
% #P4 6.0375536e-06 1.5000001 2.9906015e-08 0.020020285 2.5 2.5 2.5 2.5
% #P5 2.5 2.5 2.0525 0.8975 2 -0.5775 6 6
% #P6 1 0.15 -4.480343 -3.415989 0 
% #N 38
% #L Theta  H  K  L  Epoch  Seconds  ionm0  hubmon  ion2  ion1  temp  laser  filters  valves  TempC  topup  realt  livet  icr  vortot  vGaKb  vMoKa  vInKaNB  vInKa  vInKaPB  visual  p_peak  p_total  p_point  p_gam_s  p_del_s  p_xcen  p_xsig  p_ycen  p_ysig  p_image  pind  correct

