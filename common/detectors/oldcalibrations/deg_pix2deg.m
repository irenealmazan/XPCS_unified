function degpix = deg_pix2deg(PIXELS,DetCalib,ScnInfo,outdata)
%  
%% Correction - 
%%		RotTrans is Rotation and Translation 
%%		(like gamma on the zaxis or Delta on sevchex)
%		Rot is Rotation about an axis (conventional circle)

%% Flat detector, assume that the initial calibration is done at zero
%% 

% PIXELS [Y(row) X(col)

%        geoangles_i: [7x1 double]
%    geoangles_label: [7x5 char]
% 


DetCalib.colname 	= 'Delta';
DetCalib.RotType	= 'RotTrans';
DetCalib.colCEN	= round(mean(([1 3]-1)+116));
DetCalib.colDeg2	= 2.23;
DetCalib.col2Pix	= 100;
DetCalib.col_sgnD2P 	= +1;  % change here to neg if necessary

%  'Y' of raw image  (rows)
DetCalib.rowname 	= 'Nu';
DetCalib.RotType	= 'Rot';
DetCalib.rowDeg2	= 2.23;
DetCalib.row2Pix	= 100;
DetCalib.row_sgnD2P 	= +1;  % ch
scninfo
               docu: [3x77 char]
                 OM: [1x1 struct]
            geomode: 7
              HKL_i: [3x1 double]
        geoangles_i: [7x1 double]
    geoangles_label: [7x5 char]
        allangles_i: [36x1 double]
    allangles_label: '  Rho  Nu  Delta  Mu  Eta  Chi  Phi  DelRot                      …'
           Energy_i: 2.4000e+04





end

%%%%%%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function deg = calcDegrees(DetCalib,Nm,angm,delPix,angcen)
% choose which one to do
	if strcmp(DetCalib.RotType,'RotTrans')
		deg = calcRT(Nm,angm,delPix,angcen);
	else
		deg = calcR(Nm,angm,delPix);
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  degRT = calcRT(Nm,angm,delPix,angcen)
% Calculates the angles at a pixel for a spectrometer angle
% when the angle is created from a rotation and translation
% Assume equal pixels, but detector is flat and not curved
% and 'perpendicular' to scattered beam/direct beam
%
%   angcen is angle at detector center (degrees)
%   delPix is pixels from the center
%
%	Nm and angm are the calibrations acquired when center of detector
% 		was at zero angle

% This is exact calculation
	tandegRT 	=  delPix .* cosd(angcen).*tand(angm)./(cosd(angm).*Nm);
	degRT		=  atand(tandegRT);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  degR = calcR(Nm,angm,delPix,angcen)
% Calculates the angles at a pixel for a spectrometer angle
% when the angle is created from a rotation (circle) (normal)
% Assume equal pixels, but detector is flat and not curved,
% and 'perpendicular' to scattered beam/direct beam
%
%   delPix is pixels from the center
%
%	Nm and angm are the calibrations acquired when center of detector
% 		was at zero angle

% This is exact calculation
	tandegR 	= delPix.*tand(angm)./Nm;
	degR 		= atand(tandegR);
end
