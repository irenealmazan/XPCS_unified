function [filetextarray,error] = readdocu_files(mfilename)
% filedtextArray = readdocu_files(m_filename_str)
%
%	Designed to help self-document a parameter file by copying the 
%		file to an array. This filetextArray can then be saved within 
%		a cell or structure along with the rest of the reduced data
%
%	OUTPUT   
%		filetextArray	String Array, each line is a line in the file
%						In principle, it can be run via eval(filetextArray)
%						Or can be added to a structure to help self-document
%						the multitude of parameters that are set
%	INPUT
%		m_filename		Filename (string) of the m-file to be extracted 
%						as an array. Assume it was an 'm' file so it is 
%						not necessary to add the '.m' 
%
%	e.g., sdata.parameters = readdocu_files(RUNPARAMS)
%		
%	

if nargin<1;help readdocu_files;return;end

if isempty(which(mfilename)); 
	disp(['The parameter file ' mfilename '.m appears to be missing']);
	filetextarray = ['Error: the parameter file ' mfilename '.m was not found. Missing'];
else
	disp(['Self documenting the parameter file ' mfilename '.m ']);
	disp(['Using the file located at ' which(mfilename)]);
	
	FLAG=1;
	filetextarray = [];
	
	fid = fopen(which(mfilename),'r');
	
	while FLAG
		Aii = fgetl(fid);
		if ischar(Aii);
			filetextarray = char(filetextarray,Aii);
		else
			FLAG=0;
		end
	end
	fclose(fid)
end



