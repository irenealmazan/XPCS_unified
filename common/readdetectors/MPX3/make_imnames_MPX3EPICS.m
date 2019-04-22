function [NameMatrix]=make_imnames_MPX3EPICS(rootfilename,NUMPTS,Ndig,AREApath)
% this will make a matrix of tif filenames given roots
% Use addnames2matrix to add more to front or end of NameMatrix
% This just does the 'minimum'
%
% OUTPUT (structure)
%	NameMatrix.filenames
%	NameMatrix.fullfilenames
%	NameMatrix.specfilename	  (not used)
%
%						
%
% INPUT (plain)  assume 1 filename, 1 scannum, default is to read all points
%	rootfilename	
%	NUMPTS	if single number, go from 0 to NUMPTS (if vector, do those						
%
%   Assumption - tifs for each scan point N for scan dec1414_1 scan #S are  
%       if imname=[]
%
% 	rootfilename_NNNNNN.tif (if for example 000001, 000002, etc
%
% 

if nargin < 4;
	%[SPECpath,AREApath]=pathdisplay;
	AREApath = '/home/cthompson/DATA/2017/2017_04_mocvd_nitrideXPCS/areadata/medipix_3_EPICS/EPICS_next';
end

NUMdigitsPTNUM     = Ndig;    % NNNNN

eTIF = '.tif';
   
if length(NUMPTS(:))==1;   %assume it denotes the number from 0 to NUMPTS
	p_num = [1:NUMPTS]';
	Npts = NUMPTS;
else
	p_num = NUMPTS(:);
	Npts = length(p_num(:));
end
		
    strPTNUM    = numstring(p_num,NUMdigitsPTNUM);

% Make the other bits to string together	
    FILENAME 	= char(ones(Npts,1)*[rootfilename]);	
	FILESEP		= char(ones(Npts,1)*[filesep]);
	U		= char(ones(Npts,1)*['_']);
	TIF		= char(ones(Npts,1)*[eTIF]);
	
	FILENAMETIF			= [FILENAME,U,strPTNUM,TIF];
    FULLFILENAMETIF 	= addnames2matrix([AREApath,filesep],FILENAMETIF);
	
NameMatrix.filenames        = FILENAMETIF;
NameMatrix.fullfilenames    = FULLFILENAMETIF;
NameMatrix.points			= p_num;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%  Helper functions %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [strNNN] = numstring(N,Ndigits)
% Converts integer N to '000N' string with digits specified Ndigits
%	example, setnumstring(12,5) outputs  character string '00012'
	N=N(:);		% ensure column vector
	ten2Ndig = eval(['1e' int2str(Ndigits)]);
	strNNN	=	int2str(N + ten2Ndig);
	strNNN	=	strNNN(:,2:end);
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

