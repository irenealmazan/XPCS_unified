function  [HLs] = makeyline(x,COLOR,HA);
% [HLs] = makeyline(X,COLOR,HA); To draw a 'vertical' Y line(s) across 2d plot at position(s) X; (at maximum Z)
%  
% 	[HLs] = makeyline(X)    (uses current active plot axes) 
% 	[HLs] = makeyline(X,COLOR) (COLOR of line(s))  (default is black)
% 	[HLs] = makeyline(X,COLOR,HA) (specify HA (handle) of axes to draw lines)
%	[HLs] = makeyline(X,HA) (specify HA, use default color)
%
%	OUTPUT 
%		HLs  ;  handle(s) of the line(s)
%	INPUT
%		X is vector of positions , or more complicated (to carry through more stuff)
%		HA is handle of axes; defaults "HA=gca";  (currently active axes/figure)
%		COLOR is color of line (in RGB [1 0 0] is red) or 'r')  
%				default is black 'k' ([0 0 0])
%				may specify a column with a different color for each line ([1 1 1];[0 0 0],...)
%				or choose all to be the same color.
%
%  complementary - makexline, makeyline
%
%  		if X a structure, then it passes the X positions and the YLim that is used
%			X.Position are vectors of X positions
%			X.Lim are the 'YLim' (range) used for those x
%		if X.Lim does not exist, then program will figure out from (HA) graph

if ~isstruct(x)
	X = x;
	YLim = [];
elseif ~isfield(x,'Lim')
	YLim = [];
	X = x.Position;plot
else
	YLim = x.Lim;
	X = x.Position;
end
	

if nargin<2; COLOR = 'k'; HA = gca; end
if nargin<3 & isgraphics(COLOR,'axes');
	HA = COLOR; 
	COLOR = 'k';
elseif nargin<3
	HA = gca;
end

LINEWIDTH = 1;

if length(COLOR(:,1))<numel(X)
	if ischar(COLOR);
		COLOR = char(ones(numel(X),1)*COLOR);
	else
		COLOR = ones(numel(X),1)*COLOR;
	end
end

if ~isempty(YLim)
	YLIM = YLim;
else 
	YLIM = get(HA,'YLim');
end
	ZLIM = get(HA,'ZLim'); 

% if X is empty, numel(X) is 0 this will skip the steps and not plot anything
	for ii= 1:numel(X);
		% will work also in semilogy or loglog even if ylim includes zero
		HLs(ii) = line(HA,[0 0 0]+X(ii),[YLIM(1) (diff(YLIM).*1e-2)+YLIM(1) YLIM(2)],[0 0 0]+ZLIM(2));
		HLs(ii).LineWidth = LINEWIDTH;
		HLs(ii).Color = COLOR(ii,:);

	end
% But if X empty, will mess up output
if isempty(X); HLs=[];end

