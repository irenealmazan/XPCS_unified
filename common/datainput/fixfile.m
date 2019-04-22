function PRESET = fixfile(RUNPARAMS)
% function PRESET = fixfile(RUNPARAMS)
% can be used as >> fixfile
%  it will erase all the *.ndx files and *.cndx
%  as fixfile(RUNPARAMS) which has the ENDDATE it will not do anything past ENDDATE
%  this was originally because it also erased the data file (back when we were scp-ing them over
%  and used absence of the data file to bring over a new copy.

% if no filename is given, it will erase all the *.ndx files and *.cndx
% if filename given, it will erase its *.ndx files and
% that particular file in the data directory
% FOR EACH RUN - need to change inside this file
%       DATApath, and FINALDATE 
%

% example of dates in RUNPARAMS
%STARTDATE = '05-Dec-2017';ENDDATE ='12-Dec-2017';
 
if nargin==1; 
	eval(RUNPARAMS);
	if ~exist('ENDDATE');
		disp('we expect the run parameters file to define ENDDATE, the end of run');
		FLAG = 0;
		ENDDATE='1995-Jan-01';  % before any of our runs
	end
	% check if we are after the run date - if so, do not do anything for fixfile
	if datenum(ENDDATE)-datenum(date) <1
		FLAG = 0;
		disp('fixfile will not do anything as we are after the specified run date')
		disp('fixfile normally will erase the readspec indexing files')
		disp('and is typically using DURING the run when the data files are being updated')
		disp('and we need to reindex frequently. It is expected after ENDDATE that')
		disp('we can use the last made index files')
		disp('please see fixfile for more details if you are using it as a function');
	else
		FLAG = 1;
	end
end

if FLAG	
%% gets the datapaths currently in use by users
	[DATApath]=pathdisplay;

	% get type of computer just in case running on PCwindows
	% SINCE WE ARE USUALLY RUNNING EXPERIMENTS ON A UNIX COMPUTER - 
	% ONLY ALLOW REMOVE *.NDX (NOT THE DATA FILE) IF ON A WINDOWS COMPUTER
	if isunix
		REMOVE = 'rm ';
	else	
		REMOVE = 'del ';
	end	


     commandstring1 = ['!',REMOVE,DATApath,filesep,'*.ndx'];
     commandstring2 = ['!',REMOVE,DATApath,filesep,'*.cndx'];

	eval(commandstring1);
	disp(commandstring1);

	eval(commandstring2);
	disp(commandstring2);

end
