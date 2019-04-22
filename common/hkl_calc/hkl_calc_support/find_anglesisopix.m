function [Qxyz]=conv_isov2qxyz(HP,det_s,ScnInfo,geoscan)
%function [ANGLES]=find_anglesisopix(HP,det_s,ScnInfo,geoscan)
% Given pixX and pixY coordinates, and the detector calibration files
% that have been run correctly through detcalib_crude (eventually rotoper)
% to get the rotation operator that can convers the 

%geoscan = get_geoangles(ScnInfo,outdata,collabels,scntype);
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
%HP.Vertices = [0 0 -11;10 0 -11;0 10 -22];
%geoscan.angles = [9.6229 0 .7 36.8272];

LABELS = ScnInfo.OM.angles_label;
LABELSpix = char(det_s.Hname,det_s.Vname,'Chi','Eta');

A2PIX = chan2col(LABELS,LABELSpix);
PIX2A = chan2col(LABELSpix,LABELS);

UB = calc_UB(ScnInfo.OM);


OPER = det_s.rot.oper;
A0 = geoscan.angles(1,A2PIX);

ANGCONV = [...
		OPER, [0 0;0 0]
		0 0 1 0
		0 0 0 1];


% Just need the first point as the Chi (hardwired) in the input will have
% already have the Chi, and the other angles do not change in the scan
% so this currently does not work for hklscans or other scans

%  each rot of HV.Vertices os [pixX pixY Chi] value (if make sure to plot in chi)
NPt = length(HP.Vertices(:,1));
    ANGLES = ( ANGCONV * ...
				[(HP.Vertices)';
				zeros(1,NPt) ]  ...     
						+ A0'*ones(1,NPt))';
	
	HKL = calc_angles2hkl(ANGLES(:,PIX2A),UB,ScnInfo.Energy_i);
	
	HP.Vertices = (UB.UB * HKL.hkl')';
	
	%		|Q_phixyz> = UB|hkl(r.l.u)>
				
	

% assumes that per scan, nu and del do not change, only chi

