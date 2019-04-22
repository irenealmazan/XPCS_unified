function indexdirecotry(filelisting)
% function indexdirectory(filelisting) will do readspecscan on scan 1 for each file
%	in filelisting. This will create an index file for each scan
%	filelisting is the name of a file that has, on each line, a different specfile
%
%
%	Listing withing filelisting will look as below. It will need path if you
%	are currently residing in a different directory
%
%		/data/subdirectory/specfilename1
%		/data/subdirectory/specfilename1
%				...
%
%  Although, if you are residing in the directory with the spec data,
% and easy thing to do is
%     >> ! ls > tdirect				(this makes a file with directory listing)
%	>> ! vi tdirect				( this allows you to edit it, usually taking out
%									the tdirect file, and any other extraneous files)
%	>> indexdirectory('tdirect')		(runs this program with directory listing)
%	>> !rm tdirect				(removes the extra directory file)

% open the file with the directory listing, and read line by line
% assume each line is a different file
 fid=fopen(filelisting);
       TEST = 1;
       
        while TEST == 1;
            line = fgetl(fid);
            	if isstr(line);
            		data = readspecscan(line,1);
            	else
            		disp('end of file reached');TEST = 0;
            	end
        end
        
  fclose(fid);
 
