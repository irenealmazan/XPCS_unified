function [HL,HText] = plot_adjust(HA,colorflag,SQUEEZE,PAPERPOSITION,LINEWIDTH);
% function to set the papersize, and the fontsizes
% [HL,HText] = plot_adjust(HA,colorflag,SQUEEZE,PAPERPOSITION,LW);
%
% This is a somewhat hardwired alterations, with changes focused on what
%	has been usually convenient for the beamtime plots 
%
% OUTPUT
%	HL  handles of the lines it found
%	HText handles of the text like objects it found (title, labels)
%			Prior to R2015, get(HL(1),'Color') outputs color per line 1
%			>R2015 adds  HL.Color (for all) and HL(1).Color  (for 1)
%
% INPUT
%	HA (handle of axes, default is HA=gca to use current active)
% 	colorflag = [1 1] (default)
%				if colorflag = 0 or [0 0] or [] will not change colors or line widthgs
%				if colorflag = 1 or [1 1], will change colors and line widths
%				if colorflag = [1 0] will change colors and not line widths
%				if colorflag = [0 1] will not change colors, will change widths 
% 	SQUEEZE = [1 1 1 1] (default) no change to position or size  (also [] 
%			squeeze axes position, SQUEEZE*position
%	LW		= linewidth (default 1) 
% 
% 			A three line title looks best with a SQUEEZE = [1 1 0.93 0.88]
%			typically first 2 numbers should be 1 as they are x(y)start
%			the second two will applyt to x and y width of plots
%			Note - this squeezes the axes (so leaves more room for neighboring text)
%
%
%
LINEWIDTHdefault = [1];
SQUEEZEdefault		 = [1 1 1 1];
PAPERPOSITIONdefault = [2 2 3.75 3.5];
colorflagdefault	 = [1 1];
HAdefault = gca;
%  PAPERPOSIITONdefault = [2 2 3.5 2.8];

if nargin<5
	LINEWIDTH = LINEWIDTHdefault;
end
if nargin<4
	PAPERPOSITION =  PAPERPOSITIONdefault;
end
if nargin<3
	SQUEEZE	= SQUEEZEdefault;
end
if nargin<2;
	colorflag = colorflagdefault;
end
if nargin<1;
	HA = HAdefault; 
end

% some additional dealing in case called with [] for 'skipped' entries
if isempty(HA);HA=HAdefault;end  
if isempty(SQUEEZE);SQUEEZE = SQUEEZEdefault;end
if isempty(colorflag); colorflag = colorflagdefault;end
if isempty(LINEWIDTH); LINEWIDTH = LINEWIDTHdefault;end



% These options occasionally get changed due to specific needs for a run's plots
FONTNAME = 'Times New Roman';
		% ' Times New Roman'   'Arial'    
		% sometimes needs to be 'times' etc, see >> listfonts 
		% to get listing of all available fonts on a particular computer
		% if the msstfonts installed, there will be 'Times New Roman' and 'Arial'
		%		Arial is closest to Helvetica (embedded in postscript printers as standard font)
		% Note - matlab does NOT embed fonts in the pdf or eps so it is best to stick
		% with something that will also be available on other computers.

FONTSIZE	=	8;	% will be default for numbers, etc with axes - also, if text added later, they will be this size
CFONTSIZE	= 	8;	% used for labels, 
TFONTSIZE	=	6;	% used for titles

LINEWIDTH = LINEWIDTHdefault;
MARKERSIZE = 10;   % currently not changing marker size

% Not much to tweak below here unless you want a different default order of colors


% Find the parent figure that goes with the axes and apply the paperposition and squeezes
HF = get(HA,'parent');
	set(HF,'paperposition',PAPERPOSITION);

	set(HA,'box','on');
	set(HA,'FontSize',FONTSIZE,'FontName',FONTNAME);   % Default font with axis
	set(HA,'Position',get(HA,'Position').*SQUEEZE);

% Find the handles of the lines and texts
	HL = get(HA,'children');
	HT = get(HA,'title');
	HX = get(HA,'xlabel');
	HY = get(HA,'ylabel');
	
	HHL = [HL;HX;HY];   % will handle title separately

% When associating with the lines - the 1st line down will be lowest handle
% but bottom on legend, so it is convenience to swap order of COLORS

	COLORORDER = flip(makecolor(length(HL)));
% or use the following for totally random color sets
% 	COLORORDER = randcolor(length(HL));   

% expand out the remaining colorflag options (0 or 1)
	if length(colorflag<2) & colorflag(1)>0
			colorflag = [1 1];
	else 
			colorflag = [0 0];
	end

	% apply the changes depending on the type of graphics 
	for ii = 1:length(HHL);

		TYPE = get(HHL(ii),'Type');

		if strcmp(TYPE,'line')
			if colorflag(1)>0
		
				set(HHL(ii),'Color',COLORORDER(ii,:)); 

			end
			
			if colorflag(2)>0
				set(HHL(ii),'LineWidth',LINEWIDTH);
			end
			
		elseif strcmp(TYPE,'text');
			set(HHL(ii),'FontName',FONTNAME,'FontSize',CFONTSIZE);
		end

	end
	
	% treat title differently as having issues
	set(HT,'FontName',FONTNAME,'FontSize',TFONTSIZE,'FontWeight','bold');

end % end of function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  HELPER FUNCTIONS %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [COLORS] = makecolor(N)
% start with a fairly standard set that does not change
%		(avoiding bright yellow [1 1 0]) 
% then any additional, randomize but with colors that are not too light

	COLORS = ones(N,3);    % would be all white (too light)
	
	blue = [0 0 1];
	magenta = [1 0 1];
	cyan = [0 1 1];
	red = [1 0 0];
	green = [0 1 0];
	purple = [0.5 0 0.5];
	teal = [0 0.5 0.5];
	navy = [0 0 0.5];
	orange = [.8 .35 0];
	gray = [.6 .6 .6];
	olive = [0.7 0.7 0];

	maroon = [0.5 0 0];
	yellow = [1 1 0]; % don't use yellow usually
	
	% our creation of a standard order of some colors
	COLORS = [blue;magenta;cyan;red;green;purple;teal;navy;orange;gray;olive];


	NColor = length(COLORS(:,1));

	if N > NColor
		for ii = (NColor+1):N
	
			COLORS(ii,:) = rand(1,3);

			while mean(COLORS(ii,:))>0.9;  % keep doing it until it is not too light colored)
				COLORS(ii,:) = rand(1,3);		
			end
		end
	end

end   % end of function


%%%%%%%%%%  use this one for totally random color sets without a standard initial set
function [COLORS] = randcolor(N)
	COLORS = rand(N,3);
	
		for ii= 1:N
		
			while mean(COLORS(ii,:))>0.9;  % keep doing it until it is not too light colored)
				COLORS(ii,:) = rand(1,3);		
			end
		end
end  % end of function

