% reuse without reextracting  - presumably load from files

% unpack (eventually we should make the functions so that we are consistent and do not
% have to unpack

% check if using the same sdata and idata information
if any(id.Uid - sd.Uid) 
	disp('The data from the spec file structure and the data in the image file structure were NOT reduced together ')
	disp('This suggests that data may have been mixed together and the reduction NOT be documented correctly')
	disp('Date and time on spec data structure is ' datestr(sdata.Uid(1)) ' and the image structure is 'datestr(idata.Uid(1))]);
	disp('Other piece if ID for spec data is datestr(sdata.Uid(1)) If the dates are the same, check the Uid of both (e.g., sdata.Uid, idata.Uid)')


IInormb		= idata.IInb;
XCOLpts 	= idata.XCol;
YROWpts 	= idata.YRow;
SPECpts		= idata.Z;
timestampX	= idata.Zt;  % always 'second'close
XCOL		= idata.XColstr;
YROW		= idata.YRowstr;
Norm		= idata.Normalization;
NORMdoc		= idata.Normalizationstr;
lastframes	= idata.Valves;





% some defaults that we are not using 
THRESHusedoc = [];
FLATFIELDFLAG	= 0;




% Create Unique ID so that we can associate that the sdata and the idata have not been mixed
%Unique_ID		= [now; rand(1)];   % using datestr(Unique_ID(1) will give date and time)

%sdata.Uid		= Unique_ID;   % should check this when reusing sdata to work with altering idata 
%sdata.ReduceDate = datestr(sdata.Uid(1));
%sdata.Titlestr	= TITLEstr1

%idata.IInb 		= IInormb;
%idata.XCol 		= XCOLpts;
%idata.YRow 		= YROWpts;
%idata.Zt		= timestampX;
%idata.Ztstr		= SCNXLABEL
%idata.XColstr	= XCOL;
%idata.YRowstr	= YROW;
%idata.Normalization 	= Norm;
%idata.Normalizationstr	= NORMdoc;
%idata.Uid		= Unique_ID;
%idata.ReduceDate = datestr(sdata.Uid(1));
%idata.otherstuff = [FLATFIELDFLAG];

