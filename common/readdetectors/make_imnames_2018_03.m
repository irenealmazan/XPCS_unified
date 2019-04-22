function [NameMatrix,sdata]=make_imnames(filename,SCNNUM,readparms)
% [NameMatrix,specdata] = make_imnames(filename,SCNNUM,readingparms)
% Makes a matrix of filenames for the images based on the spec scan
%		typically spec scan needed as per point, there may be a parameter 'p_image'
%		that holds a piece of the name. (We only had a golden few years where
%		we had the placeholder for the tif point be equal to the spec point number
%	
%		this also will accommodate making a matrix for images that are not within a scan
%			typically, using ct or prpeak, will save an image. 
%			Expect a particular directory naming structure (used at beamline 12)
%			for spec point num 5 (6th point in scan) in 18th scan from scan 'filename'
%				for illustration, assume that 'p_image' for that specpoint is 7
%				because there were several images saved during filter changes
%			[path on disk]/'filename/filename_0018/filename_0018_00007.tif'
%			The tif saved from a 'ct' (need to know p_image manually, e.g., 7
%			[path on disk]/'filename/filename_000007.tif'
% 
% this will make a matrix of tif filenames from a spec scan
% Use addnames2matrix to add more to front or end of NameMatrix
% This just does the 'minimum' root file names typicall from
%		'2018_0102_1/2018_0102_1_0044/2018_0102_1_0044_00003.tif'
%			will still need to add, e.g., /DATA/DATA/2018/... 
%			that is, different locations on different computers
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
%	SCNNUM		integer value (expectation only ONE scan)
%	readingparms (structure) with fields
%		scanflag	=  1 default (a scan)   0 is from count (then ignores SCNNUM)
%		imname		= 'pil_'  or similar base tif name
%				= [] or field not exist, unique tif per out standard based on filename
%		SPECpath  spec file path 
%		AREApath  area images path (top level) (not used here), please append later with addnames2matrix
%		p_image = []   or 'pt_num'  would be (spec point)
%			= 'p_image' or other name then uses the p_image (or other) variable in data file of spec scan)
%			= [N1;N2;N3...] integer vector, not a string if having to set numbers manually 
%				(works for old situations, and to specify number for prpeak/ct (non scan) images
%			   If using several scans, then do columns for each scan, and pad with NaN
%			  Note - these MUST be vectors that are same length as the scan
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
    readparms.scanflag	= 1;	% reading a scan, = 0 is a tif without a scan (e.g., during > ct 1
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
		Npts	= sdata.Npts;
		p_num	= [];

		for ii=1:length(Npts);	

	       		if ~ischar(readparms.p_image) % to accommodate very old situations
				%disp('1')
				p_num{ii}	= readparms.p_image(:,ii);% must be column per scan
			elseif ~strcmp(readparms.p_image,'pt_num');  % then use variable in p_image string
				%disp('2')  
				p_num{ii} = round(sdata.DATA(:,chan2col(sdata.collabels,readparms.p_image)));		
			else  %% what we wish were default, use the spec point number	
				%disp('3')
				p_num{ii} = [0:Npts(ii)-1]'; 
		        end

			strSCAN(ii,:)    = numstring(SCNNUM(ii),NUMdigitsSCNNUM);
		end
		
		p_num = cat(1,p_num{:});
		p_num(isnan(p_num))=[];   % 1st option may have NaN initially		

		strPTNUM    = numstring(p_num,NUMdigitsPTNUM);

		
	else    % SCNNUM are really the particular count point numbers, e.g., from prpeak
		% also assume that they are not spread over different scans
	
		if ~ischar(readparms.p_image) % supposed to be how numbers brought over
			p_num	= readparms.p_image(:,ii);			
			Npts	= length(SCNNUM);
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
	else  % crappy non-unique names ending up like 'pil_0000' 'pil_0001'...
		FILENAMEF = char(ones(Npts,1)*readparms.imname);
	end
	
    FILENAME 	= char(ones(Npts,1)*[filename]);	
	SEP		= char(ones(Npts,1)*[filesep]);
	U		= char(ones(Npts,1)*['_']);
	SCNNN		= char(ones(Npts,1)*[strSCAN]);
	TIF		= char(ones(Npts,1)*[eTIF]);

% String them together
	
	if FLAG~=0    % from a scan
	
		if isempty(readparms.imname)
            		FILENAMETIF	= [FILENAMEF,U,SCNNN,U,strPTNUM,TIF];
            		FULLFILENAMETIF = ...
				[FILENAME,SEP,FILENAME,U,SCNNN,SEP,FILENAMEF,U,SCNNN,U,strPTNUM,TIF];
		else
			FILENAMETIF	= [FILENAMEF,strPTNUM,TIF];
			FULLFILENAMETIF = ...
				[FILENAME,SEP,FILENAME,U,SCNNN,SEP,FILENAMEF,strPTNUM,TIF];		
		end
            
	NameMatrix.pts2img	= [[0:Npts-1]'  p_num];
	NameMatrix.inputs.scnnum = SCNNUM;
		
	else   % or just a count
	
		FILENAMETIF		= [FILENAME,U,strPTNUM,TIF];
		FULLFILENAMETIF     = [FILENAME,SEP,FILENAME,U,strPTNUM,TIF];
		NameMatrix.inputs.scnnum    = SCNNUM;
		
	end

% output	
NameMatrix.filenames        = FILENAMETIF;
NameMatrix.fullfilenames    = FULLFILENAMETIF;
NameMatrix.inputs.filename  = filename;
NameMatrix.inputs.scnnum    = SCNNUM;
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
		
		specdata.DATA		=	[specdata.DATA;DATA];
		specdata.Npts		=	[specdata.Npts;Npts];
		specdata.collabels	=	[collabels];
		specdata.SCNTYPE{ii}	=	[SCNTYPE];
		specdata.SCNDATE{ii}	=	[SCNDATE];
		specdata.fileheader{ii}	=	[fileheader];
		specdata.scnheader{ii}	=	[scnheader];
		specdata.comments{ii}	=	[comments];

	end
end
