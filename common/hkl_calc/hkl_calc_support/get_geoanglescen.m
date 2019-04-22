function ANGLES = get_geoangles(ScnInfo,outdata,collabels,scntype)
% ANGLES = get_geoangles(ScnInfo,outdata,collabels,scntype)
%		extracts the geometry angles from scan and also header
%		
%	OUTPUT (structure)
%		ANGLES.angles (Npts x NAngles) matrix  
%			(if any change, will output one row per point in scan)
%		ANGLES.angles (1 x Nangles) row vector (if none show up as changed in scan

	% if scan was an hklscan, there will be many more fields carried along with ANGLES structure that describe what was needed within conversion 
%	INPUT  (will requires output from readspecscan and readspecscan_XXX_helper type program)
%		outdata (matrix) full data from the scan in spec file
%		collabels (string row) labels line from the scan in spec file
%		ScnInfo (structure) as output from readspecscan helper files
%		Requires  
%			ScnInfo.geoangles_i  (spectrometer positions from header of scan
%			ScnInfo.geoangles_label (the names of the geometry OM angles)
%			ScnInfo.OMspecial  (if working with hklscan)
%
% Typically for detector or hkl calculations, need all the geo angles

% finds which columns (NDX) in scan have motors that are geo motors
% and which ones are matching in of the geo labels (NDXgeo)

[NPTS,COLUMNS]	= size(outdata);
[geoNUM]	 	= length(ScnInfo.geoangles_i(:,1));

% if it was not a scan where the angles changed, maybe it was hklscan
if scntype(1:3)=='hkl';

	[NDX] 	= chan2col(collabels,['H';'K';'L']);
	UB 		= calc_UB(ScnInfo.OM);
	[ANGLES] = ...
            calc_hkl2angles(outdata(:,NDX),UB,ScnInfo.Energy_i,ScnInfo.OMspecial);
		

else 
            ANGLES.angles_label = ScnInfo.OM.angles_label;
% check if motors scanned include geo motors 
	[NDX,NDXgeo] 	= chan2col(collabels,ScnInfo.geoangles_label,'quiet');
	if ~isempty(NDX)
		ANGLES.angles  = ones(NPTS,1) * ScnInfo.geoangles_i(:)';
		geoangles_in_scan = outdata(:,NDX);
		% update any changing geometry motors with their values from scan
		ANGLES.angles(:,NDXgeo) = geoangles_in_scan;
	else
		ANGLES.angles = ScnInfo.geoangles_i(:)';
	end
end





