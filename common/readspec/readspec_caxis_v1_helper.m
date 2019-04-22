function [OM,scninfo,OMspecial] = readspec_caxis_v1_helper(scnheader,fileheader,comments);
% designed to pull out orientation matrix, scan information, and any special orientation information
%		specific to particular geometry 
% 2014-Dec-13 Carol Thompson cthompson@niu.edu
% First - need scnheader, fileheader, and comments from general spec reader
% Only one scan per reading 
%	>>  [outdata,outpts,scantype,scandate,ncols,collabels,scnheader,fileheader,comments] = readspecscan(name,scan)

% Extract sigma, tau, alpha, OM, EPOCH, motorpositions, motorlabels 

% OM, OMspecial
	[OM, scninfo] = extractOM(scnheader, fileheader);
	[OMspecial] = extractspecial(scnheader,fileheader);
	

end

%%%  HELPER FUNCTIONS %%%

function [OM,scninfo] = extractOM(scnheader,fileheader)
% The following will likely work for almost all more recent spec files
% and the only change needed will be the "GeoMot"line so that it specifies
% the crystallographic geometry motors and the order in the file that matches 
% the orientation matrix
% 		For example, if  the motor names and order in the OM should be tth, th, chi, phi
% 			and if the motors show up as in line below, then GeoMot = [5 1 2 3]
%  #O0 Th  Chi  Phi  zmot  TTh  Eta  Phi  DelRot
%  #O1 DelTrans
% If you are so unlucky that the geometry motors are spread over #O0, #O1, #O2, then pretend all
% positions concatenated and use those numbers.

%	In caxis, #O0     Nu Chi Eta Delta Mu
%    and in the orientation matrix, the motor order is given as 
	GeoMot = [1:7];

%  G0 and G1 lines exist in scan headers in very old spec file (circa 2000)  but G4 does not
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
	
	fileheader
	Ondx
	Pndx
	
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
	OM.a0		= {G1(19:22),fhc./G1(31)};			
	OM.a1		= {G1(25:28),fhc./G1(32)};

	scninfo.OM = OM;

% place holder for other modes, fixed beta, azimuth, theta, etc	% values going into scan (i.e., state before any action of scan) this is 'initial' values)
	scninfo.geomode	= round(G0(1));	   % should be integer
%	scninfo.Energy_i	= fhc./G4(4);	
	scninfo.HKL_i		= QHKLi(1:3);
	scninfo.geoangles_i	= [Pos(GeoMot)];
	scninfo.geoangles_label	= col2chan(Olabels,[GeoMot]);
	scninfo.allangles_i	= [Pos];
	scninfo.allangles_label	= Olabels;

% If older spec files - some parameters per scan were not kept (I assumed they assumed they were more constant)
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
%% note - since about 2011 the header of caxis has been messed up (notice all the repeated eta)
%% so that it is not saving parameters that are in the 'motor' positions. PHF rewrote the newfile 
%% command to add an option to ask user to request book and page number when 'newfile' 
%% was requested, however, it was not written compatible with 'newfile' and the headers it writes
%% out. So when this macro was listed, it clobbered the header writing macros particularly
%% if any new 'motors' were added the incorrect newfile was in place 

%%% example of output of caxis
%#F oct2714_1
%#E 1414413521
%#D Mon Oct 27 07:38:41 2014
%#C sput_141026a  User = bessrc
%#O0       Nu       Chi       Eta     Delta          Mu  TABLEVERT  TABLEHORZ
%#O1  DelS_TP   DelS_BM      mono  MONO_CHI    fpupvr    fpdnvr  MS0Bottom     MS0In
%#O2   MS0Out    MS0Top  MirBendDown  MirBendUp  MirVertUpOut  MirVertIn     Sam_Y     Sam_Z
%#O3    Sam_X  MirVertDownOut  MirTheta    MirChi   id12gap  Slit1Bot  Slit1Top  Slit1VGap
%#O4 Slit1VCen  Slit1out   Slit1in  Slit1HGap  Slit1HCen  Slit2Top  Slit2Bot  Slit2VGap
%#O5 Slit2VCen  Slit2InB  Slit2OutB  Slit2HGap  Slit2HCen    s0hgap    s0hcen    s0vgap
%#O6   s0vcen    MS1Bot    MS1Top   MS1VGap   MS1VCen    MS1InB   MS1OutB   MS1HGap
%#O7  MS1HCEN       Eta       Eta       Eta       Eta       Eta       Eta       Eta
%#C Book 92, Page 6
%
%#S 1  ascan  eta 0.829625 0.889625  30 1
%#D Mon Oct 27 07:41:26 2014
%#T 1  (Sec)
%#G0 2 0 1 0.003188262589 0.001725418314 0.9999934289 0 0 0 1
%#G1 3.905 3.905 3.905 90 90 90 1.609010322 1.609010322 1.609010322 90 90 90 2 0 1 2 0 3 22.851875 19.385 0.66775 11.2695 0 0 22.321875 43.479625 0.67375 35.881375 0 0 0.7749026676 0.7749026676
%#G3 1.600973911 -0.1606136007 -0.0001533800463 0.1606084002 1.600933065 -0.01150997611 0.001301551559 0.01143717793 1.608969146
%#G4 0.0001853693935 0.0004461859228 0.199991721 0.7749026676 1.152885409 1.121108988 -0.0006875 2.274000415 -96.87307652 0.3372 60.42435 0.9999968125 0 0 0 0 0 -180 -180 -180 -180 0 1
%#Q 0.000185369 0.000446186 0.199992
%#P0 0.001375 0 0.859625 2.274  0.00025 160.69996 73.399938
%#P1 3.299885 2.699935 7.0983204 0.1725 -3.934911 -5.065051 -0.6999998 -1.3999483
%#P2 1.400006 0.7 0.5008 0 0.50007 -1.8e-05 2.6000075 -2.111
%#P3 0.5000625 -0.50007 5.002 1.495 14.148786 -0.6825 1.7025 1.02
%#P4 -1.1925 4.18 5.8375 10.0175 0.8275 -0.0025 1.0025 1
%#P5 0.5025 1.8375 0.1625 2 -0.8375 2.0000004 2.5524764e-08 0.50000021
%#P6 0.028658151 9.475 1.625 11.3 0 444.51783 444.51783 444.51783
%#P7 444.51783 0.859625 0.859625 0.859625 0.859625 0.859625 0.859625 0.859625
%#CA Amp 1: 1.0e-06 (0.000e+00); Amp 2: 1.0e-06 (0.000e+00); Amp 3: 2.0e-08 (2.600e-10); 
%#N 62
%#L Eta  H  K  L  Epoch  Sec  ion1  Hubmon  scan_bar  Filters  Correct  Qmag  Energy  Graze  p_peak  p_total  p_del_s  p_nu_s  p_image  Vort1  Vort2  Vort3  Vort4  Vort5  Vort6  Vort7  Vort8  Vort9  Vort10  Vort11  Vort12  Vort13  Vort14  Vort15  Vort16  Vort17  Elastic  v_icr  v_ocr  v_real  v_live  v_dead  v_sum  p_xcen  p_xsig  p_ycen  p_ysig  power  tempsp  gun_on  setpw1  fwpw1  repw1  DCV1  setpw2  fwpw2  repw2  DCV2  Gamry_V  Gamry_I  ion0  p_point
%0.829625 0.000174917 0.000341993 0.199992 167 1.0000006 471110 315525 0 2 504326.14 0.32178968 15.999971 1.1228855 327241 11910091 6687308 6400759 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 132135 0 0 0 0 55.1 750 553 0 0 0 0 0 0 0 0 0 0 392883 52615
%0.831625 0.000175613 0.00034894 0.199992 168 1.0000006 469514 313326 0 2 527397.75 0.32178968 15.999971 1.1248855 343805 11818011 6632616 6349716 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 132135 0 0 0 0 55.1 750 553 0 0 0 0 0 0 0 0 0 0 391786 55022
%%%% etc