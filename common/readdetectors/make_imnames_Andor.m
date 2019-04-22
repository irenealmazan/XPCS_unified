function [FILENAMETIF,FULLNAME]=makeimnames_Andor(NUM,fileroot)
% example
%[FILENAMETIF,FULLNAME] = makeimnames_Andor([1 2 3],'2015_1118_2_')
%  FILENAMETIF is matrix - first row is  '2015_1118_2_001'
%  FULLNANE is matrix - first row is  'IMAGEpath/2015_1118_2_001'
% 	where imagepath is actually what is path to the image
%  to change the IMAGEpath - change in stpath.m

% NUM is [1 2 3] etc, makes 001 002 003 in name
% this will make a matrix of tif filenames
%		currently hardwired - only one filename, but can do more than one scan

[DATApath,IMAGEpath]=pathdisplay;
Npts = length(NUM);	

	FILENAMErootmat 	= char(ones(Npts,1)*[fileroot]);
	IMAGEpathmat		= char(ones(Npts,1)*[IMAGEpath]);
	TIF		= char(ones(Npts,1)*['.tif']);

	SSSS = NUM(:)+1e3;
		strSSSS = int2str(SSSS);
		strSSS = strSSSS(:,2:end);  % for all those zeros!
	
	FILENAMETIF	=  [FILENAMErootmat,strSSS,TIF];
        FULLNAME	=  [IMAGEpathmat,FILENAMETIF];

