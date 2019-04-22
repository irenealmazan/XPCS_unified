function DetCalib = detectorcalibrationfile

%%% Use the ROI/epics or imageJ convention for X and Y and indexing of pixels
%%% Note - indexing starts at ZERO for imageJ or the ROI/epics.
%%%  size(X, Y) and indexing for length N, ---  [0 : (N-1)]
%%% For the matlab programs that use these remember for matrices/images
%%%			(matlab) X = columns Y=rows, and indexing starts at ONE
%%% Matlab, size(Y(rows),X(columns)), and length N ---- [1 : N]
DetCalib.docu = strvcat(...
'2016-02_Pilatus file working on structure CT 2016-2 version compatible with rotated detector and non rectangular');
%% I am using XY below to remind myself it isn't in column/row 
DetCalib.pixXYSIZE	= [487 195];   % [X Y] from imageJ 
%	Xcen  Ycen
DetCalib.pixXYCEN		= [196 97];
% relative scale	  [Xscale Yscale  XYangle]
DetCalib.pixXYSCALE	= [1 1 90];

DetCalib.ADU2photon 	= -1; 
%% ADU2photon - Calibrated efficiency of 'counts' in pixel to photons; 
%%		-1 (negative number) denotes uncertain, uncalibrated
%% 		+1  (1 count in pixel is 1 photon) 
%%		(note: Pilatus are not 1 to 1 efficient at high energies.
%% 		0.5 means 1 count in pixel per 2 photons

%%  In terms of the angles (vertical) and (horizontal) 
DetCalib.HVbeam		= [0   -0.2];  % [Hor  Vert] angle where beam would be in defined center ROI (this allows scattered beam to be used, also takes account that direct beam is at delta(Vertical angle) = -rho on sevchx
DetCalib.Hname 		= 'Nu';
DetCalib.HRotType	= 'Rot';
% Vcal 2015_0208_3 #5 nu calibration scan about beam zero 
%					  [angle	Xpix	Ypix]
DetCalib.Hcal		= [ +0.76	196		32
						-0.88	196		173];


DetCalib.Vname 		= 'Delta';
DetCalib.VRotType	= 'TransRot';
% Vcal 2015_0208_3 #4 delta calibration scan about beam zero 
%					  [at_angle	Xpix	Ypix]
DetCalib.Vcal		= [ +0.80	110		97
						-1.04	275		97];
						
					

DetCalib.ADU2photon 	= -1; 

%% ADU2photon
%%	Calibrated efficiency of 'counts' in pixel to photons; 
%%		-1 (negative number) denotes uncertain, uncalibrated
%% 		+1  (1 count in pixel is 1 photon) (note: Pilatus are not 1 to 1 efficient at high energies.
%% 		0.5 means 1 count in pixel per 2 photons


%#CR ROI  1: MinX= 196  SizeX=   1 MinY=  92  SizeY=  11  <<< (VERTICAL) delta define
%#CR ROI  3: MinX= 185  SizeX=  23 MinY=  97  SizeY=   1   <<< nu (HORIZONTAL) define
%#CR ROI  4: MinX= 145  SizeX= 103 MinY=  44  SizeY= 107
end


