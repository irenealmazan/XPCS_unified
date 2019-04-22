function ...
[outdata, vtxptinfo, bins]= readVTX(name,point,reindexflag,ndxdirectory)
% [outdata, vtxptinfo, bins]= readVTX(name,point,reindexflag,ndxdirectory)
%
% outdata - column matrix of the spectrum, one scan point per column
% vtxptinfo - column matrix of the information per point 
%	(1st row is point number requested, next rows are the other information in vtx point line.
% bins - number of bins in spectrum
%
% name	- full name (with path) of the vtx file
% point - vector of points ([0 3 5] will output spectrum 0 3 and 5
% reindexflag - (optional) +1 force redo the index file
% ndxdirectory - directory (without '/' if '\' at end) for the index files
%		note - it is NOT parsed relative to the data file
%
% 25-Oct-2014 Carol Thompson   cthompson@niu.edu
%	based on readspecscan, (might be overkill for small scans, but
%		likely useful if the fluorescence detector was on for long scan.

if nargin<4; 
	ndxdirectory = [];
	filename=name;
else 
	ndxdirectory=[ndxdirectory,filesep];
	[filename,prepath] = stripfilename(name);
end
if nargin<3; reindexflag=0;end
if nargin<2; point = -1;end;   %read all points
disp(['Reading file: ' name]);

% Note: Uses mkidvtx to read entire file at the beginning, index, 
% and create index file. If already exists, each scan is located
% from the index (fseek) and read by lines (fgetl)
	
namendx = [ndxdirectory,filename '.vtxndx'];

	if 2==exist(namendx) & reindexflag~=1
  		disp(['Using existing index file ' namendx]);
  		fidi = fopen(namendx,'r');
  		fluorVTXid = fscanf(fidi,'%g %g %g %g %g %g %g',[7 inf]);
		% Get scan indexes (byte offsets)
 		idpoint = fluorVTXid(1,:);
		% Get scan numbers
  		numpoint = fluorVTXid(2,:);
		% The others are related to some info about deadtime or something not documented
  		fclose(fidi);
	else
		if 2~=exist(name) 
			disp(['Cannot find data file at  ' name])
 			outdata = -1; % use for error catching
			return
		else  %
 			[idpoint,numpoint] = mkidvtx(name,namendx);	
		end
	end

% Find index of specified spectrum point - if point requested is negative, output the spectrum at all points

if point<0;point=numpoint;end

% Make the information matrix of the points (holds what the vtx file documented about the point)
%		vtxptinfo = fluorVTXid(

outdata=[];
vtxptinfo=[];

for ii = point

idpos = find(numpoint == ii);

 if isempty(idpos)
		disp(['Fluor spectra from scan point ' num2str(point) ' not found in ' name]);
		disp(['Possibly the file has been updated so index file not up to date?']);
		disp(['    if so, readVTX with reindexflag = 1']);
    		outdata = [];
    		bins = 0;
		vtxptINFO=[];

else
% collect the point info, remember numbering starts at 0 (but matlab starts at 1)
		vtxptINFO=fluorVTXid([2:end],ii+1);

% Open up file and go to the point

  	fid = fopen(name,'r');
  	disp(['Reading spectrum data of point: ' sprintf('%d ',ii)]);
 
  	fseek(fid,idpoint(idpos(length(idpos))),'bof');
  	
	line = fgetl(fid);
	%  disp(line)
  	blanks = find(line==' ');

  	line = fgetl(fid);
	%  disp(line)
    
% Read scan data  (number of data points is actual)
    
  	data = zeros(0,1);
  	npts = 0;
  		while(1)
			line = fgetl(fid);
			if ~isstr(line),break,end;
			if length(line)==0,break,end;
			if strcmp(line(1),'P'),break,end;
			data = [data;sscanf(line,'%g',1)'];
			npts = npts + 1;
		end
	outdata = [outdata, data];

end   
%  	disp(['Length of spectra (bins): ' num2str(npts)]);
    
vtxptinfo = [vtxptinfo, vtxptINFO];
bins = npts;

 fclose(fid);

end

end  % of readspec scan main

%%%%% Helper functions %%%%%%%%%
function [idscns,numscns] = mkidvtx(name,namendx)
% [idscns,numscns] = mkidspec(name,namendx)
% used to make index of spec data file and write to disk
% name = spec data file name (string)
% 13-Aug-96 Brian Stephenson
% stephenson@anl.gov 
% elaborated from loadspec.m by Sean Brennan & Anneli Munkholm SSRL
%
% 23-Oct-2014 Carol Thompson
% cthompson@niu.edu
% Repurposed - to make index of points in the vrtx ascii file and the point number
% one spectrum per point associated with point in scan, the point line has some info
% also added some to add the rest of the info in the line
% added way for index to be elsewhere than file

disp(['Creating index file or reindexing file ' namendx]);

fid=fopen(name,'r');
chdat=fread(fid,inf,'uchar');
fclose(fid);

chdat = setstr(chdat)';

% idscns is vector of character indexes to every scan in file

idscns = findstr(chdat,'Point ');

% numscns is vector of scan numbers associated with each scan in file

numscns = zeros(size(idscns));
pointinfoi = zeros(length(idscns),2);
pointinfof = zeros(length(idscns),3);
	for ii = 1:length(idscns)
 		 numscns(ii) = sscanf(chdat(5+idscns(ii):11+idscns(ii)),'%d',1);
		 pointinfoi(ii,:) = sscanf(chdat(13+idscns(ii):30+idscns(ii)),'%i',2);
		 pointinfof(ii,:) = sscanf(chdat(31+idscns(ii):57+idscns(ii)),'%f');
	end



% Write index to disk

fidi = fopen(namendx,'w');
	if fidi < 0
 		  display('Failed to write index file; no write permission?');
	else
  		 fprintf(fidi,'%i %i %i %i %f %f %f\n',[idscns; numscns; pointinfoi'; pointinfof']);
  	 fclose(fidi);
	end

end % of mkidspec function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [filename,prepath] = stripfilename(filenamewithpath);
	SEP = strfind(filenamewithpath,filesep);
	if isempty(SEP); 
		prepath	=[];
		filename=filenamewithpath;
	else
		filename= filenamewithpath([(SEP(end)+1):end]);
		prepath	= filenamewithpath([1:max(SEP)-1]);
	end

end

%%%%%%%%%%  vtx.1 type file (ascii) looks like
%Point      0:        3        3   10.010    9.962    9.993
%0
%5
%6
%... variable number of points - same per each point in particular file, but different files may be different
%Point      1:        3        3   10.003    9.956   10.099
%


