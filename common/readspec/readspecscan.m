function ...
[outdata, outpts, scantype, scandate, ncols, collabels,fileheader,scnheader,comments]= ...
		readspecscan(name,scan)
% version 2015-07
% [outdata, outpts, scantype, scandate, ncols, collabels,fileheader,scnheader,comments]= readspecscan(name,scan)
% read in spec data file and extract specified columns from specified scan (only 1 per call)
%
% name = file name (string); scan = scalar; 
% outdata = matrix of scan data
% outpts = number of points
% scantype (string) = spec scan type
% scandate (string) = scan date and time
% ncols = number of columns in scan
% collabels (string) = labels of columns in file
% scnheader (string matrix) = full header of the scan
% fileheader (string matrix) = full header of the file (with blank lines extracted)
% comments.scn (string matrix) comments (#D, #C) before and after particular requested scan 
% comments.file (string matrix) all comments and scans commands in order (#D, #C, #S lines) 
%
% 09-Apr-98 Brian Stephenson   stephenson@anl.gov 
%
% 26-Jul-2013 Carol Thompson   cthompson@niu.edu
%	converted readspecscan to stand alone, embedding mkidspec as subfunction
%	added scnheader and fileheader outputs on request,, still backwards compatible  to gbs
%	added comments extraction on request
%
% 17-Jul-2015 Carol Thompson  cthompson@niu.edu
%	Newer fileheaders have blank lines that are not at end of header. Now accommodated.
% 

disp(['Reading file: ' name]);

% Note: Uses mkidspec to read entire file at the beginning, index, 
% and create index file. If already exists, each scan is located
% from the index (fseek) and read by lines (fgetl)
	
namendx = [name '.ndx'];

	if 2==exist(namendx)
  		disp(['Using existing index file ' namendx]);
  		fidi = fopen(namendx,'r');
  		specid = fscanf(fidi,'%g %g',[2 inf]);
		% Get scan indexes (byte offsets)
 		idscns = specid(1,:);
		% Get scan numbers
  		numscns = specid(2,:);
  		fclose(fidi);
	else
		if 2~=exist(name) 
			disp(['Cannot find data file at  ' name])
 			outdata = -1; % use for error catching
			return
		else  %
 			[idscns,numscns] = mkidspec(name);	
		end
	end

% Find index of specified scan

idpos = find(numscns == scan);

 if isempty(idpos)
		disp(['Scan ' num2str(scan) ' not found in ' name]);
    		outdata = [];
    		outpts = 0;
    		scantype = '';
    		scandate = '';
    		ncols = 0;
    		collabels = ''
 		scnheader = '';
		fileheader = '';
		comments = '';
else

% Read scan header and number of data columns

  	fid = fopen(name,'r');
  	disp(['Reading scan: ' sprintf('%d ',scan)]);
 
  	fseek(fid,idscns(idpos(length(idpos))),'bof');
  	
	line = fgetl(fid);

  	blanks = find(line==' ');

  	scantype = '';	
	if length(blanks)>=3
  		scantype = line(1+blanks(3):length(line));
		disp(['Scan type: ' scantype]);
 	 end

  	line = fgetl(fid);
 
 	scandate = '';
 	if length(line)>=3
  		if strcmp(line(1:2),'#D') 
			scandate = line(4:length(line));
			disp(['Scan date: ' scandate]);
		end
	end
 
  	ncols = 0;
  	for ii = 1:50
		line = fgetl(fid);
		if strcmp(line(1:2),'#N')
			ncols = sscanf(line(3:length(line)),'%d',1);
			disp(['Number of data columns: ' num2str(ncols)]);
			break;
		end
	end

	line = fgetl(fid);
	collabels = '';
	if length(line)>=3
		collabels = line(3:length(line));
		disp(['Column labels: ' collabels]);
  	end
    
% Read scan data  (number of data points is actual)
    
  	data = zeros(0,ncols);
  	npts = 0;
  		while(1)
			line = fgetl(fid);
			if ~isstr(line),break,end;
			if length(line)==0,break,end;
			if strcmp(line(1),'#'),break,end;
			data = [data;sscanf(line,'%g',ncols)'];
			npts = npts + 1;
		end
  	disp(['Number of data points: ' num2str(npts)]);
    
outdata = data;
outpts = npts;

	% Request read of headers for later extraction of additional metadata 
	if nargout>6

		% Specific scan header, go back one character to include the #
		fseek(fid,idscns(idpos(length(idpos)))-1,'bof');
		scnheader = []; endheader = 0;
		for ii=1:50
			line = fgetl(fid);
			if length(line)>=2 & line(1)=='#' & endheader==0
				scnheader = strvcat(scnheader,line);
			else
				endheader = 1;
			end
		end

		% File header
 	 	frewind(fid);
 		fileheader = []; endheader =0;
		for ii = 1:60
			line = fgetl(fid);
			if length(line)>=2  & line(1)=='#' & endheader==0
				if line(2)~='S';
					fileheader = strvcat(fileheader,line);
				else
					endheader = 1;
				end
			end
		end

  	end

	% request read and indexing of comments
	if nargout> 8

	namecndx = [name '.cndx']; 

	if 2==exist(namecndx)
  		disp(['Using existing comment index file ' namecndx]);
  		fidi = fopen(namecndx,'r');
  		specid = fscanf(fidi,'%g %g %g',[3 inf]);
		% Get comments indexes (byte offsets)
 		idcoms = specid(1,:);
		% Comments are before scan num in column 2 and after scan num in column 3
		% (scans command line is intercalated with byteoffset scannum scannum identical in both columns
  		numscns_pre 	= specid(2,:);
		numscns_post    = specid(3,:);
  		fclose(fidi);
	else  
		[idcoms,numscns_pre,numscns_post] = mkidspecC(name,idscns,numscns);		
	end

	% Specific scan header
		comments.scn = '';
		comments.file = '';
		
		if ~isempty(idcoms)
			for ii=1:length(idcoms);
				fseek(fid,(idcoms(ii)-1),'bof');
				line = fgetl(fid);
					comments.file = strvcat(comments.file,line);
				if (numscns_pre(ii)==scan | numscns_post(ii)==scan)
					comments.scn = strvcat(comments.scn,line);
				end
			end
		end

	end

 fclose(fid);

end

end  % of readspec scan main

%%%%% Helper functions %%%%%%%%%

function [idscns,numscns] = mkidspec(name)
% [idscns,numscns] = mkidspec(name)
% used to make index of spec data file and write to disk
% name = spec data file name (string)
% 13-Aug-96 Brian Stephenson
% stephenson@anl.gov 
% elaborated from loadspec.m by Sean Brennan & Anneli Munkholm SSRL

namendx = [name '.ndx']; 

disp(['Creating index file ' namendx]);

fid=fopen(name,'r');
chdat=fread(fid,inf,'uchar');
fclose(fid);

chdat = setstr(chdat)';

% idscns is vector of character indexes to every scan in file

idscns = findstr(chdat,'#S ');

% numscns is vector of scan numbers associated with each scan in file

numscns = zeros(size(idscns));
	for ii = 1:length(idscns)
 		 numscns(ii) = sscanf(chdat(2+idscns(ii):6+idscns(ii)),'%d',1);
	end

% Write index to disk

fidi = fopen(namendx,'w');
	if fidi < 0
 		  display('Failed to write index file; no write permission?');
	else
  		 fprintf(fidi,'%i %i\n',[idscns; numscns]);
  	 fclose(fidi);
	end

end % of mkidspec function

%%%%%%%%%%%%%%%%%%%%%%

function [idcoms,numscns_pre,numscns_post] = mkidspecC(name,idscns,numscns)
% [idcom] = mkidspecC(name)
% pick out comments, (and any date remarks)
% output filestructure is  [charbyte_of_comment    beforethisscan   afterthisscan]

namecndx = [name '.cndx']; 
Nscn = length(numscns);

disp(['Creating comment index file ' namecndx]);

fid=fopen(name,'r');
chdat=fread(fid,inf,'uchar');
fclose(fid);

chdat = setstr(chdat)';
lastcharacter = length(chdat(:));

% idcoms will be vector of character indexes to every comment and date statement in file

idcom  = findstr(chdat,'#C');
iddates = findstr(chdat,'#D');

idcoms 		= sort([idcom iddates]');
numscns_pre 	= zeros(size(idcoms));
numscns_post	= zeros(size(idcoms));


		for ii=1:Nscn
			idscnii	= idscns(find(numscns == (ii)));

			idposiin 	= find(numscns == (ii-1));
				if isempty(idposiin)
					idscniin = idcoms(1)-1;
				else
					idscniin = idscns(idposiin);
				end
			idposiip	= find(numscns == (ii+1));
				if isempty(idposiip)
					idscniip=max([idcoms(end) idscns(end)]) +1;
				else 
					idscniip = idscns(idposiip);
				end

			numscns_pre	= ([idcoms > idscniin] & [idcoms < idscnii]) .*ii  + numscns_pre;
			numscns_post 	= ([idcoms < idscniip] & [idcoms > idscnii]).*ii  + numscns_post;
		end

idcmat = [	idcoms(:)  	numscns_pre(:)	numscns_post(:) 
		idscns(:)		numscns(:)	numscns(:) ];

% intercalate so scan line is in between relevant comments
idcmat = sortrows(idcmat,1);

idcoms 		= idcmat(:,1);
numscns_pre	= idcmat(:,2);
numscns_post	= idcmat(:,3); 

% Write comments index to disk
fidi = fopen(namecndx,'w');

	if fidi < 0
 		  display('Failed to write comments index file; no write permission?');
	else
  		 fprintf(fidi,'%i %i %i\n',[idcmat']);
  	 fclose(fidi);
	end 

end % of mkidspecC function

%%%%%%%% All scans in Basic Spec (Certified Scientific) program scan have at least  following in scan header
%%  sample file header
%% number of #O lines (motor names) should  match number of #P lines (position of motor before scan starts)
%
% #F nov2299_3
% #E 943320429   %Epoch -
% #D Mon Nov 22 19:27:09 1999
% #C FeCrNi  User = mocvd
% #O0    Delta     Theta        Mu     Gamma   GamTran    GamRot        sh        sv
%
%%    sample scan header
%% number of #G and #P lines and each entry parsing will depend on diffractometer  and spec version running.
%
% #S 14  tseries 2000 0.5
% #D Sun Jul  9 20:00:50 2000
% #T 0.5  (Seconds)
% #G0 3 0 0 0 0 1 0 0 0 0 0 0
% #G1 3.905 3.905 3.905 90 90 90 1.60901 1.60901 1.60901 90 90 90 1 0 0 0 1 0 60 30 0 0 0 0 60 30 0 -90 0 0 1.54 1.54
% #G2 0
% #Q 0.000251922 -0.0251895 -1.0021
% #P0 15.5 9.19 90 -0.573 -2.6666667 -1.6715265 -2.5142857 0.62870206    % there may be #P1, #P2 positions of motors
% #N 8
% #L Time  Epoch  Seconds  filters  ic2  ssd2  ic3   ic3     %labels 

