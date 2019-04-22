function [LABELS,INTEREST] = load_scans_nplot(filenames,sSCNNUMS,RUNPARAMS,DOCU,FORCE);
eval(RUNPARAMS);  % eventually uses monname, detname, VALVE, the names
SCNNUMS = eval(sSCNNUMS);
% note - this one needs readspechelper (see helper file pointer in RUNPARAMS to extract info out of headers of spec)

if nargin<5;
	FORCE.FLAG=0;FORCE.X='none';
	FORCE.FLAGdet=0;FORCE.det='none';
end

[DATApath,IMAGEpath] = pathdisplay;                                

if length(filenames(:,1))==1;
	filenamematrix = char(ones(length(SCNNUMS),1)*filenames);
else
	filenamematrix = filenames;
end

IY=[];X=[];H=[];K=[];L=[];EPOCH=[];SCNPTS=[];HL=[];clear VALVE;

for ii = 1:length(SCNNUMS);

	fullname = [DATApath,filesep,filenamematrix(ii,:)];
	[DATA, outpts, SCNTYPE, SCNDATE, ncols, LABELS,fileheader,scnheader,surroundingcomments] = ...
		readspecscan(fullname,SCNNUMS(ii));
	
	%% specific to headers from this spectrometer 
	%% (in princinple, this has pointer set up in RUNPARAMS to point to helper file for 
	%% particular spectrometer (CT has ones for caxis, and for sevchex, and for zaxis)
	[OM,scninfo,OMspecial] = readspechelp(scnheader,fileheader,surroundingcomments);
			% note    HKL at start of scan is    scninfo.HKL_i	
    
    % I use SCNFLAG here, so that I can tweak here if 
	% I want to make more complicated, and here it
	% is easier to add whether I want h or l scan if it
	% says hklscan

	if ~isempty(findstr(SCNTYPE,'tseries'))|~isempty(findstr(SCNTYPE,'loopscan'))|~isempty(findstr(SCNTYPE,'trigger'));
	
	SCNFLAG = 'GROW';
	[HLii,IYii,Xii,Hii,Kii,Lii,Xname,Ylabeltext,VALVE] = plotdata_image(DATA,SCNFLAG,LABELS,RUNPARAMS,scninfo.HKL_i,FORCE);
	
	elseif ~isempty(findstr(SCNTYPE,'hklscan'));
	
	SCNFLAG = 'HKLS';
	[HLii,IYii,Xii,Hii,Kii,Lii,Xname,Ylabeltext,VALVE] = plotdata_image(DATA,SCNFLAG,LABELS,RUNPARAMS,scninfo.HKL_i,FORCE);
	
	else	
	SCNFLAG = 'NONE';
	[HLii,IYii,Xii,Hii,Kii,Lii,Xname,Ylabeltext,VALVE] = plotdata_image(DATA,SCNFLAG,LABELS,RUNPARAMS,scninfo.HKL_i,FORCE);
	
	end
	
	SCNPTSii = length(IYii);
	VALVES(ii) = VALVE;

	HL	= [HL;HLii];
	SCNPTS  = [SCNPTS;SCNPTSii];
	X 	= [X;Xii];
	IY	= [IY;IYii];
	H	= [H;Hii];
	K	= [K;Kii];
	L	= [L;Lii];
	EPOCH	= [EPOCH; DATA(:,chan2col(LABELS,'Epoch'))];



end

% if all SCNPTS same then likely may want a grid, then put in grid instead for output
% This needs to be fixed to use if we return multiple columns of IY (say for multiple detectors)
% I think right now I am not using this for grids
if 0
	if all(diff(SCNPTS)==0)
		IY = reshape(IY,SCNPTS(1),length(SCNPTS));
		H = reshape(H,SCNPTS(1),length(SCNPTS));
		K = reshape(K,SCNPTS(1),length(SCNPTS));
		L = reshape(L,SCNPTS(1),length(SCNPTS));
		EPOCH = reshape(EPOCH,SCNPTS(1),length(SCNPTS));
		X = reshape(X,SCNPTS(1),length(SCNPTS));
	end

	% in case people put scnnums in column
	SCNNUMS = reshape(SCNNUMS,1,length(SCNNUMS));

end

	TITLESTRING = strvcat(...
		[pfilename(filenames(1,:)),' #',(sSCNNUMS),' at ',SCNDATE],...
		[pfilename(SCNTYPE)],...
		[pfilename(DOCU)]);
	title(TITLESTRING)

	% the following adjusts paper size of the plot, and
	% adjusts other prettiness, so that the three lines
	% in the title will appear on the plot
	%plot_adjust(HA,colorflag,SQUEEZE,PAPERPOSITION,LINEWIDTH);
	plot_adjust(gca,[0 1],[1 1 0.95 0.95],[2 2 4.5 3.5],0.5);
	%set(gca,'Userdata',[X H K L]);

	% add the valve switches etc at end (otherwise they
	% if they are empty, nothing will show up on plot)
	%if ~isempty(valveswitch);
	%	Xswitch = X(valveswitch);
	%	makeyline(Xswitch,[0 0 1]);
	%else
	%	Xswitch = [];
	%end

% The VALVES.Position are the spec point number (so in matlab indexing, need to point to N+1)
% Our X is long column, so need to accommodate a translate for each scan
	TRANS = 0;
	for ii = 1:numel(SCNNUMS);

		VALVEStrans = VALVES(ii);
			if ~isempty(VALVES(ii).Position);
				VALVEStrans.Position = X(VALVES(ii).Position + TRANS + 1);   % turn into plot's x position
			end  % otherwise, remi
		HL_valveii = makeyline(VALVEStrans,'b');
		set(HL_valveii,'DisplayName','Valve');
		disp([filenamematrix(ii,:) ' Scan #',int2str(SCNNUMS(ii)),': Spec point number at valve switches: ' num2str(VALVES(ii).Position(:)') ' : ' Xname,' at ' num2str(VALVEStrans.Position')]);
		TRANS = TRANS + SCNPTS(ii);
	end

 

INTEREST.IY =IY;
INTEREST.X = X;
INTEREST.H = H;
INTEREST.K = K;
INTEREST.L = L;
INTEREST.EPOCH = EPOCH;
INTEREST.SCNPTS = SCNPTS;
INTEREST.HL = HL;
INTEREST.title = TITLESTRING;
INTEREST.Xname = Xname;
INTEREST.LABELS = Ylabeltext;

end

%%%%%%%%%%%%  HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%
function [HL1,IY,X,H,K,L,Xname,LABELTEXT,VALVESWITCH] = plotdata_image(DATA,SCNFLAGS,LABELS,RUNPARAMS,HKLinit,FORCE)

% 2017-08 CT converted findstr (deprecated) to strfind
%% 2017-10 CT note - need to go through this - clean up is XLABELTEXT and Xname redundant in use later?
% 2017-10 CT added the valves (used to have these a long long time ago!
% 2018-03 CT added FORCE.FLAGdet = -1 (for filter corrected but counting to monitor)



% defaults for detectors
eval(RUNPARAMS);
monC 	= chan2col(LABELS,monname);
secC	= chan2col(LABELS,timename);

if nargin<6, 
	FORCE.FLAG=0;FORCE.X='none';
	FORCE.FLAGdet=0;FORCE.det=detname;
end

if FORCE.FLAGdet;   %anything nonzero
	detname = FORCE.det;   %otherwise, detname default from RUNPARAMS
else  % if 0, then use the default detector
	FORCE.det = detname;
end
	
IYuncorrected = DATA(:,chan2col(LABELS,detname));

	%		= 1 (det*filters/mon)  = 2 (det/mon)  = 3 (raw)  = 4 
	% I have taken out the moncpsave but want to get it back in since that helps compare between filters
	
	if FORCE.FLAGdet==1
	filtC	= chan2col(LABELS,filtername);
	CORRECTION = filter_correction(DATA(:,filtC),filter1,filtercorr)./DATA(:,monC);
	YLABELTEXT = ['Int = ' pfilename(col2row(FORCE.det)) '(filtercorrected )/(' pfilename(monname),')'];
	
	elseif FORCE.FLAGdet==2
	CORRECTION = 1.0 ./DATA(:,monC);
	YLABELTEXT = ['Int = ' pfilename(col2row(FORCE.det)) '/ ' pfilename(monname)];	
	
	elseif FORCE.FLAGdet==3
	CORRECTION = 1;
	YLABELTEXT = ['Int = ' pfilename(col2row(FORCE.det))];

	else  % default current like =2

	CORRECTION = 1.0 ./DATA(:,monC);
	YLABELTEXT = ['Int = ' pfilename(col2row(FORCE.det)) '/ ' pfilename(monname)];	

	end

% May have several columns if several detectors selected
IY = IYuncorrected .* (CORRECTION * ones(1,length(IYuncorrected(1,:))));

MaxI		= max(IY);
MinI		= min(IY);
if length(IY)>15;
    StartI		= mean(IY(1:5));
else 
    StartI  = mean(IY);
end

if 	~isempty(strfind(SCNFLAGS,'GROW'));

		etimeC	= chan2col(LABELS,etimename);
		X 	= DATA(:,etimeC); 
		H = HKLinit(1)*ones(size(X));
		K = HKLinit(2)*ones(size(X));
		L = HKLinit(3)*ones(size(X));
		HKLdocu = sprintf(': H=%7.4f K=%7.4f L=%7.4f',HKLinit(1),HKLinit(2),HKLinit(3));
		if ~FORCE.FLAG
			Xname = etimename;
			XLABELTEXT = [etimename,' [sec] at ',HKLdocu];
		else 
			if strfind(FORCE.X,'specpoint')
				X = [0:length(IYuncorrected)-1]';
				XLABELTEXT = ['spec point num '];
				Xname = XLABELTEXT;
			else
			Xname = FORCE.X;
			X = DATA(:,chan2col(LABELS,Xname));
			XLABELTEXT = [FORCE.X];
			end
		end
        


elseif 	~isempty(strfind(SCNFLAGS,'HKLS'));
		HC = chan2col(LABELS,hname);
		KC = chan2col(LABELS,kname);
		LC = chan2col(LABELS,lname);
		H = DATA(:,HC);
		K = DATA(:,KC);
		L = DATA(:,LC);
		HKL = [H K L];
		HKLdocu = [	'h'
				'k'
				'l'];
		if ~FORCE.FLAG | isempty(chan2col([hname '  ' kname '  ' lname],FORCE.X)) 
		%not forced to be a particular HKL 
			% here I am trying to make a way to figure out
			% which was most 'interesting' hkl from which
			% one varied the most
			DIFFH = max(H) - min(H);
			DIFFK = max(K) - min(K);
			DIFFL = max(L) - min(L);
			% can use SCNHKL to pick out correct xlabel
			% and the correct HKL column to plot

			[dum,SCNHKL] = max([DIFFH  DIFFK  DIFFL]);
			SCNHKL = SCNHKL(1);
			X = HKL(:,SCNHKL);
			Xname = [HKLdocu(SCNHKL)];
		else
			if strfind(FORCE.X,'specpoint')
				X = [0:length(IYuncorrected)-1]';
				Xname = ['spec point num '];
			else
			Xname = FORCE.X;
			X = HKL(:,chan2col([hname '  ' kname '  ' lname],Xname));
			end
		end
		XLABELTEXT = Xname;

else		% probably an ascan of some type

		if FORCE.FLAG 
			% allow option to plot in spec point if FORCE.X='specpoint';
			if strfind(FORCE.X,'specpoint')
				X = [0:length(IYuncorrected)-1]';
				Xname = ['spec point num '];
			else
				XC = chan2col(LABELS,FORCE.X);
				Xname = col2chan(LABELS,XC);
				X = DATA(:,XC);
			end
		else XC=1;
			% the following should work if in all other scans
			% the columns following the X will be H K L
			% that seems to be the case for most of the ascans
			% this documents all the moving motors
			%Xname = LABELS(1:(findstr(LABELS,hname)-1));
			Xname = col2chan(LABELS,XC);
			X = DATA(:,XC);
		end
		
		H = DATA(:,chan2col(LABELS,hname));
		K = DATA(:,chan2col(LABELS,kname));
		L = DATA(:,chan2col(LABELS,lname));
		HKLdocu = sprintf(': [H=%7.4f K=%7.4f L=%7.4f]',H(1),K(1),L(1));
	
		XLABELTEXT = [Xname ' at ',HKLdocu];
		
end


% PLOTTING SEQUENCE
% separated out to make it easier to change and tweak
% or for later split of this program into data analysis and
% data plotting

if 	~isempty(strfind(SCNFLAGS,'GROW'));

	HL1 = line(X,IY);
	makecolor(HL1);


elseif 	~isempty(strfind(SCNFLAGS,'HKLS'));

	HL1 = line(X,IY);
	set(gca,'Yscale','log');
	makecolor(HL1);

else

	HL1 = line(X,IY);
	makecolor(HL1);

end

%	makecolor(HL1);
	ylabel(YLABELTEXT);
	xlabel(XLABELTEXT);

LABELTEXT.X = XLABELTEXT;
LABELTEXT.Y = YLABELTEXT;


XLIMl = [min(X) max(X)];	 TESTx = diff(XLIMl);
	if TESTx~=0; XLIMl = XLIMl.*[0.98 1.02];end
YLIMl = [min(IY(:)) max(IY(:))];  TESTy=diff(YLIMl);  
    if TESTy~=0; YLIMl = YLIMl.*[0.98 1.02];else YLIMl = YLIMl + [-1 1];end
set(gca,'Xlim',XLIMl,'Ylim',YLIMl);

% for those runs that have some valves that we are tracking 2017-1005
% 2018 08  pulling this out so that the valves are plotted at the end

if exist('valvename')==1
	valveC  = chan2col(LABELS,valvename);
	if ~isempty(valveC);
%	valveswitch = find(diff(DATA(:,valveC))~=0); %this would be matlab index before valve switch works for mocvd
	VALVESWITCH.Position = find(abs(diff(DATA(:,valveC)))>=diffvalve); 
	VALVESWITCH.Lim = YLIMl;
	%as calculated 
	%	this valveswitch index is equivalent to the index of the spec point number at valve switch
	%	or if indexing in matlab, the valve switch would be one point later
	%		so using to point to an index in a matlab arrays, use as valveswitch+1
	%		but if using to point to the spec data point, use as valveswitch
	else
	VALVESWITCH.Position = [];
	VALVESWITCH.Lim = YLIMl;
	end
end

% If there are valve switches defined, put a line where they are
% haven't checked for multiple plots yet (i.e. multiple scan plots CT 2017-1005
% not implemented yet, but should get ready for more than one valve being 
%if ~isempty(valveswitch);
%	Xswitch = X(valveswitch);
%	makeyline(Xswitch,[0 0 1]);
%else
%	Xswitch = [];
%end

end
%%%%%%%%%%  HELPER FUNCTIONS %%%%%%%%%%%%%
function [COLORS] = makecolor(HL)

COLORS = [1 1 1];
for ii=1:length(HL)
	while mean(COLORS)>0.8;  % keep doing it until it is not too light colored)
	COLORS = rand(1,3);
	end

	set(HL(ii),'color',COLORS);
	COLORS = [1 1 1];  % reset if necessary to go through loop again
end

% note could do all at once with HL(:).Color = [r g b;r g b;...];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RowLabels] = col2row(Labels)
RowLabels=['['];
[NR,NC] = size(Labels);
	if NC ~=1
		for ii=1:NR;
			Labelii = Labels(ii,find(Labels(ii,:)~=':'));
			RowLabels = [RowLabels,' ',Labelii];
		end
		RowLabels = [RowLabels,']'];
	else
			Rowlabels = Labels;
	end
end



