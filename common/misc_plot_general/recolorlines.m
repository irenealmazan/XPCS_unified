function [HL,COLORS] = recolorlines(HA,avoidDisplayNames)
% Find lines in plot and randomize to a new set of colors trying to not be too light
% [HL,COLORS] = recolorlines(HA,avoidDisplayNames);  
%		if HA (handle of axes) not given, uses current active axes/plot HA=gca;
% OUTPUT  
%		HL  handle(s) of N 'line' we found associated with axes
%		COLORS  [Nx3] matrix [R G B] colors used for each line
% INPUT
%		HA (optional) if not given it will use currently active axes/plot
%				(i.e., default is HA=gca;)
%		avoidDisplayNames (string array) although better is cell array  (optional)
%			don't change lines with HL.DisplayName  (I use it for valves)
%
%		note - recolorlines will not 'find' the lines associated with the legend
%			but they will automaticall recolor to match  the plot lines
%		for runs, unless specified, it will not recolor 'Valve'

if nargin<1;HA=gca;end

if nargin <2 ; 
	aDN = [];
	aDN = 'Valve';  % for purposes of not having to worry about it, if really want to color valves,
					% then just put some other string as placeholder in nargin
elseif  ~iscell(avoidDisplayNames) & ~isempty(avoidDisplayNames)
	% probably given a matrix of names that needs to be turned to cells and trimmed for easy strcmp
	aDN = strtrim(num2cell(avoidDisplayNames,2));
end

HL = get(HA,'children');

	for ii = 1:length(HL);

		TYPE = get(HL(ii),'type');

		if strcmp(TYPE,'line');
		
			if [~isempty(aDN) &  any(strcmp(aDN,get(HL(ii),'DisplayName')))];
				% something matched, so do not change
				COLORS(ii,:) = get(HL(ii),'Color');  % will keep color
			else % go ahead and change color
				COLORS(ii,:)=randcolor(HL(ii));
			end
		end

	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [COLORS] = randcolor(HL)

COLORS = [1 1 1];

while mean(COLORS)>0.9;  % keep doing it until it is not too light colored)
	COLORS = rand(1,3);
end

	set(HL,'color',COLORS);

end

