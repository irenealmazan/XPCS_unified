function [HL,COLORORDER] = showrois(ROIS,HF,LW,COLORORDER)
% [HL,COLORORDER]= showrois(ROIS,HF,LW,COLORSROIS)  HF figure, LW linewidth, COLOR roi lines
% as written Matlab 2014b or earlier will still work as well as later versions
% OUTPUT
%	HL(N). 	graphics handles HL(1) for ROI1 of the ROI lines
%				Note ROIs that are single points will plot as a circle
%	COLORORDER	the colororder used
% INPUT
%	ROIS	2D array of N x 4- ROIS should be in [Xmin Xmax Ymin Ymax;...] 
%				in whatever units the plot is in that is, this one works simply from axes
%	HF		Figure handle to use (will use current figure, assumes only one axis
%				(figure handle 
%					probably can accept legacy (e.g., integer N to denote figure N)
%	LW		Line width (default 2)
%	COLORORDER	(N,3) RGB  [R G B;... values 0 to 1 ([1 0 0] is red)
%			default color order is set up (but can change if tweak inside)
%
% CT 2018-08 added ability to see 'single' point ROIS
%


	if nargin<2;	
		HF=gcf;LW=2;writeROIflag=1;COLORORDER=[];
	elseif nargin<3
		LW=2;writeROIflag=0;COLORORDER=[];
	elseif nargin<4
		writeROIflag=0;COLORORDER=[];
	else
		writeROIflag=0;
	end

	ndxrois = [1:length(ROIS(:,1))];

	if isempty(COLORORDER);

		[COLORORDER] = makecolororder(length(ROIS(:,1)));

	end

	Xl = ROIS(:,1);
	Xh = ROIS(:,2);
	Yl = ROIS(:,3);
	Yh = ROIS(:,4);

%	Instead of making axis active by specifying figure(HF)
%		could probably do  HA = get(HF,'children') not sure which version
%		of matlab the following line would work
%	then HL = line(HA,([Xl Xh Xh Xl Xl]'),([Yl Yl Yh Yh Yl]'));
	figure(HF);
	
	HL = line(([Xl Xh Xh Xl Xl]'),([Yl Yl Yh Yh Yl]'));

%	find the single point ROIS	
	DOTS = find(~[[diff(ROIS(:,[1 2])')] + [diff(ROIS(:,[3 4])')]]);
	
for ii=ndxrois;

	set(HL(ii),'LineWidth',LW,'Color',COLORORDER(ii,:),'LineStyle','--');
	set(HL(ii),'DisplayName',['ROI #',int2str(ii)]);
		if any(ii==DOTS);
			set(HL(ii),'Marker','o');
		end
	% If no one is using R2014a or earlier, can start using following protocols
	% HL(ii).LineWidth	=	LW
	% HL(ii).Color		= COLORORDER(ii,:);
end

if writeROIflag;
STRING = strvcat('ROIs shown are the following',...
' ROIindex : minX maxX  minY maxY  ',...
' ',...
['ROInum : MYROIS = [...'],...
int2str([ndxrois' ROIS]),'  ];');

disp(STRING)
end

end   % of function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Helper Functions%%%%%%%%%

function [COLORORDER] = makecolororder(N)
	magenta = [1 0 1];
	cyan = [0 1 1];
	yellow = [1 1 0]; % don't use yellow usually in line plots, but over the surface images might work
	red = [1 0 0];
	blue = [0 0 1];  % often the images are blue background 
	green = [0 1 0];
	gray = [.6 .6 .6];
	maroon = [0.5 0 0];
	olive = [0.7 0.7 0];
	teal = [0 0.5 0.5];
	navy = [0 0 0.5];
	purple = [0.5 0 0.5];
	orange = [.8 .35 0];
	
	COLORORDER = [magenta;cyan;green;red;purple;teal;navy;orange;gray;olive;blue];

	NColor = length(COLORORDER(:,1));

% if more than the colors above, the later ROI colors will be random (and change plot to plot)
	if N > NColor
		for ii = (NColor+1):N
	
			COLORS(ii,:) = rand(1,3);

			while mean(COLORS(ii,:))>0.9;  % keep doing it until it is not too light colored)
				COLORS(ii,:) = rand(1,3);		
			end
		end
	end
end
