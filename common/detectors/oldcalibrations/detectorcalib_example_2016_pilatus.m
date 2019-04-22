function DetCalib = detectorcalibrationfile

%DetCalib.docu
%DetCalib.SIZE = [195 487];  % [Y=rows X=columns] = size(image)  
%DetCalib.ADU2photon = [-1] %uncalb

%  'X' of raw image  (columns) what it is when seeing without any flips etc from imageJ
%DetCalib.colname 	= 'Delta';
%DetCalib.coldir 	= '+z'; % at 0,0,0,0,0, axis of spectrometer associated with Pos direction of colname 
%DetCalib.colCEN	= round(mean(([1 3]-1)+116));
%DetCalib.colDeg2	= 2.23;
%DetCalib.col2Pix	= 100;
%DetCalib.col_sgnD2P 	= +1;  % change here to neg if necessary

%  'Y' of raw image  (rows)
%DetCalib.rowname 	= 'Nu';
%DetCalib.rowdir 	= '+x';  % at 0,0,0,0,0, axis of spectrometer associated with Pos direction of rowname 
%DetCalib.rowCEN	= round(mean(([1 7]-1)+100)); 
%DetCalib.rowDeg2	= 2.23;
%DetCalib.row2Pix	= 100;
%DetCalib.row_sgnD2P 	= +1;  % change here to neg if necessary

#CR ROI  1: MinX= 192  SizeX=   3 MinY=  91  SizeY=  11

DetCalib.docu = strvcat(...
'2016-0113 CT test calibration file working on structure');

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
DetCalib.coldir 	= '+z'; % at 0,0,0,0,0, axis of spectrometer associated with direction of positive Angle 
%%			associated with colname. Note - detector pixels still may be go negative with respect to this.
DetCalib.colCEN	= round(mean(([1 3]-1)+192));
DetCalib.colDeg2	= 2.23;
DetCalib.col2Pix	= 100;
DetCalib.col_sgnD2P 	= +1;  % change here to neg if necessary

%  'Y' of raw image  (rows)
DetCalib.rowname 	= 'Nu';
DetCalib.RotType	= 'Rot';
DetCalib.rowdir 	= '+x';    % at 0,0,0,0,0, what axis of spectrometer is direction, positive delta moves detector in +x 
DetCalib.rowCEN	= round(mean(([1 11]-1)+91)); 
DetCalib.rowDeg2	= 2.23;
DetCalib.row2Pix	= 100;
DetCalib.row_sgnD2P 	= +1;  % change here to neg if necessary

%% ADU2photon
%%	Calibrated efficiency of 'counts' in pixel to photons; 
%%		-1 (negative number) denotes uncertain, uncalibrated
%% 		+1  (1 count in pixel is 1 photon) (note: Pilatus are not 1 to 1 efficient at high energies.
%% 		0.5 means 1 count in pixel per 2 photons
DetCalib.ADU2photon 	= -1; 

%#CR ROI  1: MinX=   1  SizeX= 487 MinY=   1  SizeY= 195
%#CR ROI  2: MinX= 116  SizeX=   3 MinY= 100  SizeY=   7
%#CR ROI  3: MinX= 106  SizeX=  23 MinY=  90  SizeY=  27
%#CR ROI  4: MinX=  76  SizeX=  83 MinY=  60  SizeY=  87
end


