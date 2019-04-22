function detcalib = GaN_2013_07
%	detcalib.docu = various documentations about choices
%	detcalib.rowang2pix = deg/pixels for e.g., gam (use NEG if pixel number gets smaller with bigger angle, + opposite);
%	detcalib.colang2pix =  deg/pixels for e.e., del
%		these will typically be the same unless the pixel size is rectangular)
%		for tifs, X is column, Y is row
%	detcalib.columnL = 'Delta'  for example, the spec detector 'motion' name that maps onto the pixels
%	detcalib.rowL = 'Gamma'
%		to calibrate  this - move detector arm a fixed amount (i.e., beam fixed)
%			if beam spot moves to lower pixel with increasing angle of detector - then angle/pix is POSITIVE
%			that is, it means that higher pixel number on detector is at a higher angle (if beam were moving)
%			at 0,0,0,0,0  detector   with respect to spectrometer  +Z (up) (gam +), +X (inboard) (del +), 
%					not on detector, but other axis is +Y (downstream) (right hand axis)
%			at 0,0,0,0,0 detector  '+x'  in matrix (columns) is typically  +Z (w.r.t. spectrometer), '+y' in matrix (rows) is +X
%					
%	detcalib.ADU2photon = detector 'values' in change to counts  
%			not necessarily well calibrated number, use 1 if unknown
%	detcalib.rowcenpix = center (gam) in image - pixel where is 'center' (i.e, if at 0,0
%				where would beam hit)
%	detcalib.colcenpix =  image center (e.g., gam)
%	detcalib.SIZE  = [rows(Delta)  columns(Gamma)]    % note  rows 'y'(Delta), columns 'x' (Gamma) in detector direction
%										in spectrometer direction at 0,0,0,0,... [X,Z], 

detcalib.docu = strvcat(...
'2013-0329 CT nitride mocvd run 2013-07  - (no great deg2pix determination) much confusion with arms',...
'the active area is near the bottom, worked to change apertures, note images are tif 32 bit',...
'ROI 3 and 4 changed occasionally, but p\_point ROI 2 was mostly consistent [116 118 100 106] (3x7)',...
' 	note at beamline, ROI expressed, 116 (start of X) 3 (width)  100 (start of Y) 7 (width)',...
'	tif is 497(X,columns) x 195(Y,rows)   without transpose, X (columns) Gamma: Y (rows) Delta',...	
'       and from saved Tiffs, +Gamma +ang2pix, +Delta +ang2pix',...
'	that is, if there are two spots on detector, the spot at higher pixel is at higher Gamma or Delta');


%	as I use, if a beam hits at pixel A2, is it higher or lower in gam(del) compared to a beam that hits at pixel A1
%	if A2 is bigger than A1 pixel, and if beam hitting at pixel A2 is at a higher gam(del) than at A1 then I define ang2pix positive)

%Ddeg2pix  = +2/85;	% keep this positive for ease of calculating - usually the same both directions
Ddeg2pix =  +2.23/100;   % 2.23 deg per 100 pixels , actually pretty close to 2deg/85pixel from older runs

% Do this for 'raw' tif - what are rows and columnes does it read out when doing imread
% we often later have to transpose to 'understand' it (e.g., often gamma in real space is 'up'
% and delta is sideways)
detcalib.SIZE = [195 487];  % [rows columns] of matrix in matlab that the tif or whatever outputs
% that is, what does >> size(image) output for raw image when read by matlab (no transpose, no nothing)
% Note - when looking at image, X is columns, Y is rows [Nrow,Ncol]=size(image) [NY,NX]=size(image) as seen in plotting

%  'X' of raw image  (columns)
detcalib.colname 	= 'Gamma';
detcalib.colCEN	= round(mean(([1 3]-1)+116));
detcalib.coldeg2pix 	= +Ddeg2pix;  % change here to neg if necessary
detcalib.coldir 	= '+z'; % at 0,0,0,0,0, what axis of spectrometer is direction, positive gamma is +z moves detector in +z

%  'Y' of raw image  (rows)
detcalib.rowname 	= 'Delta'; 
detcalib.rowCEN	= round(mean(([1 7]-1)+100)); 
detcalib.rowdeg2pix 	= +Ddeg2pix;  % change here to neg if necessary
detcalib.rowdir 	= '+x';    % at 0,0,0,0,0, what axis of spectrometer is direction, positive delta moves detector in +x


detcalib.ADU2photon 	= 1; 

%#CR ROI  1: MinX=   1  SizeX= 487 MinY=   1  SizeY= 195
%#CR ROI  2: MinX= 116  SizeX=   3 MinY= 100  SizeY=   7
%#CR ROI  3: MinX= 106  SizeX=  23 MinY=  90  SizeY=  27
%#CR ROI  4: MinX=  76  SizeX=  83 MinY=  60  SizeY=  87
end


