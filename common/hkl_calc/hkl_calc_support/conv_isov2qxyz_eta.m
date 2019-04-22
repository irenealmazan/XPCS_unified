function [HKL]=conv_isov2qxyz(HP,det_s,ScnInfo,geoscan)
%function [HKL,ANGLES]=conv_isov2qxyz(det_s,ScnInfo,geoscan)
%function [ANGLES]=find_anglesisopix(HP,det_s,ScnInfo,geoscan)
% Given pixX and pixY coordinates, and the detector calibration files
% that have been run correctly through detcalib_crude (eventually rotoper)
% to get the rotation operator that can convers the 

%geoscan = get_geoanglescen(ScnInfo,outdata,collabels,scntype);
%ScnInfo out of readspecscan helper file
%det_s from detectorcalibration file run through detcalib_crude or better
% HP is handle of the isosurface patch

% This is the operator to change Xpixels,Ypixel pairs to angles
% det_s.rot.oper = [D2PH 0 ; 0 D2PV] * (DetRotOper');
%% that is 
%  HORangle  =   oper   *  pixX
%  VERangle                pixY

% currently hardwired for single (chi) geo motor scans for caxis
%

%DUmmy vertices to test
%HP.Vertices = [0 0 -30.16;50 0 -30.16;0 50 -30.16;
%				0 0 -30.16;50 0 -30.16;0 50 -30.16;
%				0 0 -10.26;50 0 -10.26;0 50 -10.26];
%geoscan.angles = [9.6229 0 .7 36.8272];

LABELS = ScnInfo.OM.angles_label;
LABELSpix = char(det_s.Hname,det_s.Vname,'Eta','Rho','Mu','Chi','Phi');

%% chan2col(collabels,choice) gives the 'choice' column in collabels
A2PIX = chan2col(LABELS,LABELSpix);
PIX2A = chan2col(LABELSpix,LABELS);

UB = calc_UB(ScnInfo.OM);


OPER = det_s.rot.oper;
A0 = geoscan.angles(1,:);   


%hardwired to take out the motor already in pixz (Rho  Nu  Delta  Mu  *Eta  Chi  Phi) 
A0 = [A0(1:4) 0 [A0(6:7)]];

ANGCONV = [...
		OPER, [0 0 0 0 0;0 0 0 0 0]
		0 0 0 0 0 0 0
		0 0 0 1 0 0 0
		0 0 0 0 1 0 0
		0 0 0 0 0 1 0
		0 0 0 0 0 0 1];

Nmot = length(A0)-3;

% Just need the first point as the Chi (hardwired) in the input will have
% already have the Chi, and the other angles do not change in the scan
% so this currently does not work for hklscans or other scans

%  each rot of HV.Vertices os [pixX pixY Chi] value (if make sure to plot in chi)
NPt = length(HP.Vertices(:,1));
    ANGLES = ( ANGCONV * ...
				[(HP.Vertices)';
				zeros(Nmot,NPt) ]  ...     
						+ A0(A2PIX)'*ones(1,NPt))';
	HKL = calc_angles2hkl(ANGLES(:,PIX2A),UB,ScnInfo.Energy_i);
	
% if non cubic, need to work with this first to get orthogonoal and etc
%	HP.Vertices = (UB.UB * HKL.hkl')';

%	HP.Vertices = HKL.hkl;
		HP.Vertices = (UB.UB*HKL.hkl')';
	
	%		|Q_phixyz> = UB|hkl(r.l.u)>
				
	

% assumes that per scan, nu and del do not change, only chi

