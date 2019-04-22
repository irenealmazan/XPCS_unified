function [NameMatrix,sdata]=make_imnames_2018_03(filename,SCNNUM,readparms)
% [NameMatrix,specdata] = make_imnames_2018_03(filename,SCNNUM,readingparms)
% Changing so it will work with several scans (same root file) at once
% this will make a matrix of tif filenames from a spec scan
% Use addnames2matrix to add more to front or end of NameMatrix
% This just does the 'minimum'
%
% OUTPUT (structure)
%	NameMatrix.filenames
%	NameMatrix.fullfilenamesm
%	NameMatrix.pts2num		[specpt  tifnum]	
%	NameMatrix.inputs		carries inputs (additional fields, filename,SCNNUMS,
%	
%	specdata.{fields} all the outputs of readspecscan as fields, 
%			If not reading spec file to get num of points or p_image variable
%			then specdata.DATA is empty and no other fields exist
%						
%
% INPUT (plain)  assume 1 filename, 1 scannum, default is to read all points
%	filename	specfilename
%	SCNNUM		integer value (can be more than one scan (but later programs need to change)
%				exception, for ct/prpeaks SCNNUM is ignored, since tifs are under filename
%	readingparms (structure) with fields
%		scanflag	=  1 default (a scan)   0 is from count (then ignores SCNNUM)
%		imname		= 'pil_'  or similar base tif name
%				= [] or field not exist, unique tif per out standard based on filename
%		SPECpath  spec file path 
%		AREApath  area images path (top level) (not used here), please append later with addnames2matrix
%		p_image = []   or 'pt_num'  would be (spec point)
%			= 'p_image' or other name then uses the p_image (or other) variable in data file of spec scan)
%			= [N1; N2; N3...] integer vector, not a string if having to set numbers manually
%				this is for images created in the prpeak/ct or non scan image
%				(still need filename since these are saved in directories that depend on spec file)
%			  Also, for very old situations, we have to specify MANUALLY the 'pt' numbers
%			  (this was when the controllers were sequentially assigning numbers to the tifs, but 
%			 the numbers simply advanced throughout the file_) (But still saved within filename structures)
%			  In this old situation, if using several scans, then each scan per columns
%			 each column for each scan, and since columns must be same length, padded with NaN

%		ending	= 'tif' (default) will add this ending (e.g.,) '.tif' to end
						
%
%   Assumption - tifs for each scan point N for scan dec1414_1 scan #S are  
%       if imname=[]
%
% 	{AREApath}/dec1414_1/dec1414_1_SSSS/dec1414_1_SSSS_NNNNN.tif
%
%       if imname= 'pil_' (or similar depending)
%
% 	{AREApath}/dec1414_1/dec1414_1_SSSS/pil_SSSS_NNNNN.tif  << for older crappy nontraceable names
% 	{AREApath}/dec1414_1/dec1414_1_SSSS/pil_NNNNN.tif  << FOR caxis crappy non-traceable names/non unique
%
%   where the SSSS (scan number) and NNNNN (point number) are padded with zeros
%
%   Assumption - tifs for 'ct' or prpeak (where SCNNUM are numbers held in p_image
%
%       {AREApath}/dec1414_1/dec1414_1_NNNNN.tif  << this is for caxis and others
% 
% 

		
% Note - when things are stable, they followed the above assumptions,
% When the pilatus or area detector was first installed, or upgraded
% then there would be a run or several where things would be not consistent
% 
% Note, to construct namematrix there may be blanks at end of some lines
%   (if necessary when using, can use deblank later when using the strings
%           deblank(FILENAMETIF(1,:)); for examples
%


NUMdigitsSCNNUM    = 4;    % SSSS
NUMdigitsPTNUM     = 5;    % NNNNN

% Default is to assume we are setting up to read image set with scan
if nargin<3;
    readparms.scanflag	= 1;  % =1 comes from a scan, =0 for a 'ct' or similar 
    [SPECpath,AREApath]	= pathdisplay;   % use defaults paths already set up
    readparms.SPECpath	= SPECpath;
    readparms.AREApath	= AREApath;
    readparms.imname	= [];	% default for zaxis and sevchex, caxis 'pil_'
    readparms.p_image	= 'p_image';  %=[]  just use the point number
    readparms.ending	= '.tif';
    
end


FLAG    = readparms.scanflag;  
eTIF	= readparms.ending;

if isempty(readparms.p_image);
	readparms.p_image = 'pt_num';
end
	
NPTS=0;
	
	if FLAG~=0	%from a scan

		sdata	= DATAstruct(filename,SCNNUM,readparms);
		Npts	= sdata.Npts;   % Number of points in a scan 
		p_num	= [];

		
		for ii=1:length(Npts);	% {cell} for each scan number which may have different 

	       		if ~ischar(readparms.p_image) % to accommodate very old situations (not tested)
				
				p_num{ii}	= readparms.p_image(:,ii);% must be column per scan

			elseif ~strcmp(readparms.p_image,'pt_num');  % then use variable in p_image string

				% Need to change chan2col to accept cells eventually to be able to combine scans with different labels
				sdata.collabels = cat(1,sdata.collabels{1});   % this cuts it to 1st and makes it non  cell
				p_num{ii} =  round(sdata.DATA(:,chan2col(sdata.collabels,readparms.p_image)));	
	
			else  %% what we wish were default, use the spec point number	
				
				p_num{ii} = [0:Npts(ii)-1]'; 
		        end

			strSCAN(ii,:)    = numstring(SCNNUM(ii),NUMdigitsSCNNUM);
		end
		p_num = cat(1,p_num{:});
		p_num(isnan(p_num))=[];   % 1st option may have NaN padding initially		

		strPTNUM    = numstring(p_num,NUMdigitsPTNUM);

		
	else    % SCNNUM are really the particular count point numbers, e.g., from prpeak
		% also assume that they are not spread over different scans
		% haven't really tested this
	
		if ~ischar(readparms.p_image) % supposed to be how numbers brought over
			p_num	= readparms.p_image(ii);	% only one filename, so only one column or row		
			Npts	= length(p_num);
			%strSCAN	= 'ct';
			strSCAN		= [];
			strPTNUM = numstring(p_num,NUMdigitsPTNUM);

		end
	
		sdata.DATA=[];
		sdata.filename	= filename;
		sdata.comments	= 'These are images from cts or prpeaks; monitor or times are unknown';
		sdata.Npts	= length(p_num);
	end

% Make the bits to string together

	if isempty(readparms.imname)  % unique filenames associated with files, etc
		FILENAMEF = char(ones(Npts,1)*[filename]);
	else  % crappy non-unique or tracable names
		FILENAMEF = char(ones(Npts,1)*readparms.imname);
	end
	
    	FILENAME 	= char(ones(Npts,1)*[filename]);	
	SEP		= char(ones(Npts,1)*[filesep]);
	U		= char(ones(Npts,1)*['_']);
	SCNNN		= char(ones(Npts,1)*[strSCAN]);
	TIF		= char(ones(Npts,1)*[eTIF]);

% String them together to make the rows in the name matrix
	
	if FLAG~=0    % from a scan
	
		if isempty(readparms.imname)
            		FILENAMETIF	= [FILENAMEF,U,SCNNN,U,strPTNUM,TIF];
%					   dec1414_1_SSSS_NNNNN.tif 
            		FULLFILENAMETIF = ...
				[FILENAME,SEP,FILENAME,U,SCNNN,SEP,FILENAMEF,U,SCNNN,U,strPTNUM,TIF];
%				 dec1414_1/dec1414_1_SSSS/dec1414_1_SSSS_NNNNN.tif
		else
			FILENAMETIF	= [FILENAMEF,strPTNUM,TIF];
%					   pil_NNNNN.tif
			FULLFILENAMETIF = ...
				[FILENAME,SEP,FILENAME,U,SCNNN,SEP,FILENAMEF,strPTNUM,TIF];
%				 dec1414_1/dec1414_1_SSSS/pil_NNNNN.tif		
		end
            
	NameMatrix.pts2img	= [[0:Npts-1]'  p_num];
	NameMatrix.inputs.scnnum = SCNNUM;
		
	else   % or just a count
	
		FILENAMETIF		= [FILENAME,U,strPTNUM,TIF];
%				 	   dec1414_1_NNNNN.tif
		FULLFILENAMETIF     = [FILENAME,SEP,FILENAME,U,strPTNUM,TIF];
%				 	dec1414_1/dec1414_1_NNNNN.tif
		NameMatrix.inputs.scnnum    = p_num;
		
	end

% output	
NameMatrix.filenames        = FILENAMETIF;    
NameMatrix.fullfilenames    = FULLFILENAMETIF;
NameMatrix.inputs.filename  = filename;  % at present, while we can have several scan numbers, all need to be in one filename
%NameMatrix.inputs.scnnum    = SCNNUM;
NameMatrix.inputs.readparms = readparms;

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
function [specdata] = DATAstruct(filename,SCNNUM,readparms)

	% assume same type of scan (so same collabels) 
	% and same run, so same specpath for all files

	SPECp   = readparms.SPECpath;
	specdata.SCNNUM		=	SCNNUM;
	specdata.filename	=	filename;

	specdata.DATA		=	[];
	specdata.Npts		=	[];
	specdata.SCNTYPE	=	[];
	specdata.SCNDATE	=	[];
	specdata.collabels	=	[];    
	specdata.fileheader	=	[];
	specdata.scnheader	=	[];
	specdata.comments	=	[];

	for ii=1:length(SCNNUM);		
		SPECp		= readparms.SPECpath;
		[DATA,Npts,SCNTYPE,SCNDATE,ncols,collabels,fileheader,scnheader,comments] ...
				= readspecscan([SPECp,filesep,filename(ii,:)],SCNNUM(ii));

		specdata.Npts		=	[specdata.Npts;Npts];
		specdata.DATA		=	[specdata.DATA;DATA];  % wanted to turn to cells but need to change more
		specdata.collabels{ii}	=	[collabels];  % scans have to have same type of data
		specdata.SCNTYPE{ii}	=	[SCNTYPE];
		specdata.SCNDATE{ii}	=	[SCNDATE];
		specdata.fileheader{ii}	=	[fileheader];
		specdata.scnheader{ii}	=	[scnheader];
		specdata.comments{ii}	=	[comments];
		% I'd like to change to 'cell' the collabels to accomodate if different scans had slightly
		% different stuff (sometimes things change quickly) but that requires change of
		% chan2col and col2chan (might make them easier if I require collabels to be cells
		% and also in this program, would need to specify which 'scan' as go through identifying columns for p_image

	end
end
