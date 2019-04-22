function  [HLs] = makexline(y,COLOR,HA);
% [HLs] = makexline(Y,COLOR,HA); To draw a 'horizontal' X line(s) across 2d plot at position(s) Y (at Zmax)
%  
% 	[HLs] = makexline(Y)    (uses current active plot axes) 
% 	[HLs] = makexline(Y,COLOR) (COLOR of line(s)) e.g., [0 0 0;1 1 1;...] (default is black [0 0 0])
% 	[HLs] = makexline(Y,COLOR,HA) (specify HA (handle) of axes to draw lines)
%	[HLs] = makeyline(Y,HA) (specify HA, use default color)
%
%	OUTPUT 
%		HLs  ;  handle(s) of the line(s)
%	INPUT
%		Y is vector of positions , 
%		HA is handle of axes; defaults "HA=gca";  (currently active axes/figure)
%		COLOR is color of line (in RGB [1 0 0] is red) or 'r')  
%				default is black 'k' ([0 0 0])
%				may specify a column with a different color for each line ([1 1 1];[0 0 0],...)
%				or choose all to be the same color.
%
%  complementary - makexline, makeyline 

if ~isstruct(y)
	Y = y;
	XLim = [];
	disp('g')
elseif ~isfield(y,'Lim')
	XLim = [];
	Y = y.Position;
	disp('h');
else
	XLim = y.Lim;
	Y = y.Position;
end

if nargin<3; HA = gca;end
if nargin<2; COLOR = 'k'; HA = gca; end
if ishandle(COLOR); HA = COLOR; COLOR = 'k';end

LINEWIDTH = 1;

if length(COLOR(:,1))<length(Y)
	if ischar(COLOR);
		COLOR = char(ones(length(Y),1)*COLOR);
	else
		COLOR = ones(length(Y),1)*COLOR;
	end
end

if ~isempty(XLim)
	XLIM = XLim;
else 
	XLIM = get(HA,'XLim');
end
	ZLIM = get(HA,'ZLim');

for ii= 1:numel(Y);
	% will work also in semilogx or loglog even if xlim includes zero
	HLs(ii) = line(HA,[XLIM(1) (diff(XLIM).*1e-2)+XLIM(1) XLIM(2)],[0 0 0]+Y(ii),[0 0 0]+ZLIM(2));
	HLs(ii).LineWidth = LINEWIDTH;
	HLs(ii).Color = COLOR(ii,:);

end
