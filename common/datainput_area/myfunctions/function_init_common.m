
%	function_init_common - these may change whenever we play with new stuff even if we have the data readin
	% not sure if these are used later, but in case
	Nr = length(YROWpts);
	Nc = length(XCOLpts);
	Nz = length(SPECpts);

% these flags are used, and may change 
if ~LOGFLAG
	DOCUInt = ['Inorm',THRESHusedoc];
else
	DOCUInt = ['log10(Inorm',THRESHusedoc,')'];
end

% This FLATFIELDFLAG should not change unless we re-reduce the IInormb (images) idata
if FLATFIELDFLAG
	DOCUInt = [DOCUInt,'[FF]'];
end

TITLEstr1 = char(...
		[pfilename(HDFfilename),' #', SCNstr,' : ', sdata.SCNDATE],...
		[pfilename(sdata.SCNTYPE)],...
		[pfilename(DOCU0),' ',pfilename(DOCUscan)]);

TITLEstrshort = char(...
		[pfilename(HDFfilename),' #',SCNstr,' : ',DOCUInt,' [',pfilename(DOCU0),']'],...
		[pfilename(sdata.SCNTYPE), ' : ',pfilename(DOCUscan)]);
		
	% I expect that I will make readspechelp files compatible with cells for the above
	[OM,scninfo,OMspecial]	= readspechelp(sdata.scnheader,sdata.fileheader,sdata.comments);
	INFOndx 				= chan2col(scninfo.allangles_label,INFOinterest);
	INFOpositions 			= scninfo.allangles_i(INFOndx);
		DUM=[];
		% put into a line
		for ii=1:length(INFOpositions);
			DUM = [DUM,num2str_round(INFOpositions(ii)) ' : '];
		end
	INFOstr = char(INFOinterest,DUM);   % currently added onto 0 and 5
	% INFOstr = [];    % turn this off quickly if desired 
	
	
% decide what we are plotting for 'z' axis (time or some other variable)z

	if FORCE.FLAG; %Xstepcol=1;end   % use scan variable	 
		if strfind(FORCE.X,'specpoint')
			SCNXLABEL = 'spec point # ';
			Xstepcol=-1;
		else
			SCNXLABEL	= FORCE.X;
			Xstepcol = chan2col(sdata.collabels,FORCE.X); 
				if isempty(Xstepcol);disp('ERROR: asking for Xaxis not present in file');end
		end	
	else
		Xstepcol = 1;
		SCNXLABEL = col2chan(sdata.collabels,1); 
	end
	
	if strncmp(SCNXLABEL,etimename,4); %usually the elapsed time name is 'Time'
		SCNXLABEL	= [SCNXLABEL,' [\Delta sec] from SpecFile '];  % in pixirad, 
		Xsteps	= timestampX-timestampX(1);
	elseif Xstepcol<1;
		Xsteps	= [1:length(timestampX)]-ImageJ;     
	else
		Xsteps	= sdata.DATA(:,Xstepcol);
	end
