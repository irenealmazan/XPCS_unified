function DetCalib = detectorcalibrationfile

%DetCalib.docu
%DetCalib.SIZE = [195 487];  % [Y=rows X=columns] = size(image)  
%DetCalib.ADU2photon = [-1] %uncalb

%  'X' of raw image  (columns) what it is when seeing without any flips etc from imageJ
%DetCalib.colname 	= 'Delta';
%DetCalib.coldir 	= '+z'; % at 0,0,0,0,0, axis of spectrometer associated with Pos direction of colname 
%DetCalib.colCEN	= round(mean(([1 3]-1)+116))-1;
%DetCalib.colDeg2	= 2.23;
%DetCalib.col2Pix	= 100;
%DetCalib.col_sgnD2P 	= +1;  % change here to neg if necessary

%  'Y' of raw image  (rows)
%DetCalib.rowname 	= 'Nu';
%DetCalib.rowdir 	= '+x';  % at 0,0,0,0,0, axis of spectrometer associated with Pos direction of rowname 
%DetCalib.rowCEN	= round(mean(([1 7]-1)+100))-1; 
%DetCalib.rowDeg2	= 2.23;
%DetCalib.row2Pix	= 100;
%DetCalib.row_sgnD2P 	= +1;  % change here to neg if necessary

% For center, take the MEDM ROI values
%   START   Ns;    SIZE  NSIZE, 
%   ROI enumerating all pixels (in Pilatus/ImageJ indexing) 
%		   [1:NSIZE]-1 + Ns   
%	ROI enumerating all pixels (in matlab matrix indexing 
%		   ( [1:NSIZE]-1 + Ns ) + 1

%% Pixels in epics (ROI) /imageJ indexing starts at [0,0] 
%% Pixels in Matlab  indexing starts at [1,1]

DetCalib.docu = strvcat(...
'2016-02_Pilatus file working on structure');
% #CR ROI  1: MinX= 192  SizeX=   3 MinY=  91  SizeY=  11

%	as I use, if a beam hits at pixel A2, is it higher or lower in gam(del) compared to a beam that hits at pixel A1
%	if A2 is bigger than A1 pixel, and if beam hitting at pixel A2 is at a higher gam(del) than at A1 then I define ang2pix positive)

% Do this for 'raw' tif - what are rows and columns does it read out when doing imread
DetCalib.SIZE = [195 487];  % [Y=rows X=columns] = size(image) 

%% Please - here you must consider the RAW image (and how RAW data image is read by matlab or imageJ)
%
%% Beware - Do not use the 'imageJ' pixel readings from the image opened (epics, start_viewer)
%%		on the data control screen. Typically  there are a lot of idiosyncratic remappings (rotated, mirrored, etc)
%%		
%% Open a image file directly from a different version of imageJ  or with matlab.
%% 
%% The epics ROI's are set up using the RAW image. It will be a real mistake if anyone changes that
%% since then it is impossible (since no one ever documents the remappings) for a user later to look at
%% the RAW data and figure out which pixels were actually where the ROI was located.

%% Fix beam, move detector - if POSITIVE detangle Beam Spot moves MORE POSITIVE pixel, then sgnD2P negative
%%			sgnD2P is related to thinking about fixed detector, if beam moves)

%%
%%  Image is a matrix; 
%%		in matrix the top left hand corner is [1,1]
%%		the bottom right hand corner is [NY,NZ]   (Y (up down)is number of rows, X (right left) is number of columns)
%%		In the example below, size is [3 4], 
%%			pixel A is indexed [1,1], pixel L is indexed [3,4]  
%%   
%%	[A B C D
%%   E F G H
%%   I J K L]
%%

% Do this for 'raw' tif - what are rows and columns does it read out when doing imread
DetCalib.SIZE = [195 487];  % [Y=rows X=columns] = size(image)  

%  'X' of raw image  (columns) what it is when seeing without any flips etc from imageJ
DetCalib.colname 	= 'Delta';
DetCalib.RotType	= 'RotTrans';
DetCalib.coldir 	= '-z'; % at 0,0,0,0,0, axis of spectrometer associated with direction of positive Angle 
%%			associated with colname. Note - detector pixels still may be go negative with respect to this.
DetCalib.colCEN	= round(mean(([1 1]-1)+196))+1;
DetCalib.colDeg2	= 1.7;
DetCalib.col2Pix	= 146;
DetCalib.col_sgnD2P 	= +1;  % change here to neg if necessary

%  'Y' of raw image  (rows)
DetCalib.rowname 	= 'Nu';
DetCalib.RotType	= 'Rot';
DetCalib.rowdir 	= '+x';    % at 0,0,0,0,0, what axis of spectrometer is direction, positive delta moves detector in +x 
DetCalib.rowCEN	= round(mean(([1 1]-1)+97))+1; 
DetCalib.rowDeg2	= 1.7;
DetCalib.row2Pix	= 146;
DetCalib.row_sgnD2P 	= +1;  % change here to neg if necessary

%% ADU2photon
%%	Calibrated efficiency of 'counts' in pixel to photons; 
%%		-1 (negative number) denotes uncertain, uncalibrated
%% 		+1  (1 count in pixel is 1 photon) (note: Pilatus are not 1 to 1 efficient at high energies.
%% 		0.5 means 1 count in pixel per 2 photons
DetCalib.ADU2photon 	= -1; 

%#CR ROI  1: MinX= 196  SizeX=   1 MinY=  92  SizeY=  11  <<< delta define
%#CR ROI  2: MinX= 185  SizeX=  23 MinY=  84  SizeY=  27
%#CR ROI  3: MinX= 185  SizeX=  23 MinY=  97  SizeY=   1   <<< nu define
%#CR ROI  4: MinX= 145  SizeX= 103 MinY=  44  SizeY= 107
end


