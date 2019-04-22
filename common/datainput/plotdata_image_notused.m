function [HL1,IY,X,H,K,L,Xname,LABELTEXT] = plotdata(DATA,SCNFLAGS,LABELS,RUNPARAMS,HKLinit,FORCE)
PRESET = 1;
% 2017-08 CT converted findstr (deprecated) to strfind
%% 2017-10 CT note - need to go through this - clean up is XLABELTEXT and Xname redundant in use later?
% 2017-10 CT added the valves (used to have these a long long time ago!
% 2018-03 CT added FORCE.FLAGdet = -1 (for filter corrected but counting to monitor)



% defaults for detectors
eval(RUNPARAMS);
monC 	= chan2col(LABELS,monname);
secC	= chan2col(LABELS,timename);

% for those runs that have some valves that we are tracking 2017-1005
valveswitch=[];
if exist('valvename')==1
	valveC  = chan2col(LABELS,valvename);
	if ~isempty(valveC);
	valveswitch = find(abs(diff(DATA(:,valveC)))>=diffvalve); %this would be spec point number before valve switch
	end
end

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

	
	if FORCE.FLAGdet==1
	filtC	= chan2col(LABELS,filtername);
	CORRECTION = filter_correction(DATA(:,filtC),filter1,filtercorr)./(DATA(:,secC).*DATA(:,monC));
	YLABELTEXT = ['Int = ' pfilename(col2row(FORCE.det)) '(filtercorrected )/(time*' pfilename(monname),')'];
	
	elseif FORCE.FLAGdet==-1
	filtC	= chan2col(LABELS,filtername);
	CORRECTION = filter_correction(DATA(:,filtC),filter1,filtercorr)./DATA(:,monC);
	YLABELTEXT = ['Int = ' pfilename(col2row(FORCE.det)) '(filtercorrected )/(' pfilename(monname),')'];
	
	elseif FORCE.FLAGdet==2
	CORRECTION = 1.0 ./DATA(:,monC);
	YLABELTEXT = ['Int = ' pfilename(col2row(FORCE.det)) '/ ' pfilename(monname)];	
	
	elseif FORCE.FLAGdet==3
	CORRECTION = 1;
	YLABELTEXT = ['Int = ' pfilename(col2row(FORCE.det))];

	else
	CORRECTION = 1.0 ./(DATA(:,secC).*DATA(:,monC));	
	YLABELTEXT = ['Int = ' pfilename(col2row(FORCE.det)) '/(time*' pfilename(monname),')'];
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


XLIM = [min(X) max(X)];	 TESTx = diff(XLIM);
	if TESTx~=0; XLIM = XLIM.*[0.98 1.02];end
YLIM = [min(IY(:)) max(IY(:))];  TESTy=diff(YLIM);  
    if TESTy~=0; YLIM = YLIM.*[0.98 1.02];else YLIM = YLIM + [-1 1];end
set(gca,'Xlim',XLIM,'Ylim',YLIM);

% If there are valve switches defined, put a line where they are
% haven't checked for multiple plots yet (i.e. multiple scan plots CT 2017-1005
if ~isempty(valveswitch);
	Xswitch = X(valveswitch);
	makeyline(Xswitch,'b');
else
	Xswitch = [];
end

end  % function

%%%%%%%%%%  HELPER FUNCTIONS %%%%%%%%%%%%%
function [COLORS] = makecolor(HL)

COLORS = [1 1 1];
for ii=1:length(HL)
	while mean(COLORS)>0.8;  % make sure not too light colored)
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
