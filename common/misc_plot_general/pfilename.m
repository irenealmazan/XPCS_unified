function PRETTYFILEmat = pfilename(SPECFILEmat)
% Inserts a "\" in front a all "_"  in input string (making it suitable for titles/labels of plots) 
% PRETTYname = pfilename(filename)
%
% OUTPUT   
%		string (or string matrix) with '\_' in place of all '_'
% INPUT
%		string or a string matrix  of filenames
%
% The latex interpreter assumes '_' is command for subscript
%	however, '\_' will appear as '_' in titles and labels of plots

PRETTYFILEmat = [];
for jj = 1:length(SPECFILEmat(:,1));

    % Find the offending characters
    SPECFILE = SPECFILEmat(jj,:);
    IDX = find(SPECFILE == '_');
    PRETTYFILE = [];
    % Insert \_ for  all _  found
    for ii = 1:length(SPECFILE)
	    if any(ii==IDX)
		    PRETTYFILE = [PRETTYFILE,'\',SPECFILE(ii)];
	    else
		    PRETTYFILE = [PRETTYFILE,SPECFILE(ii)];
	    end
    end
PRETTYFILEmat = strvcat(PRETTYFILEmat,PRETTYFILE);
end
