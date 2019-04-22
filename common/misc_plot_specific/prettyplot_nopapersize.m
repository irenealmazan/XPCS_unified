function [PRESET] = plotstuff(HF);
% applied to figure and works if more than one figures (subplot)
% etc for the current plot
% call by
%	plotstuff(HF)
%	where HA is the handle for the current figure
%	HA = gca;
% This one has a SQUEEZE parameter that was set
% up to allow a three line title 
% do all the stuff except change the papersize (useful if have bigger plot and want to 

if nargin<1;HF=gcf;end

%PAPERPOSITION = [2 2 3.5 2.8];
FONTSIZE = 10;
FONTNAME = 'times';
LINEWIDTH = 0.8;
MARKERSIZE = 10;
COLORORDER = ['b';'g';'r';'c';'k';'m'];

% need to do modulo here, but needed quicklyfor last plot of evening
% WILL BOMB OUT IF STILL MORE LINES THAN COLORS LISTED BELOW
COLORORDER = [COLORORDER;COLORORDER;COLORORDER;COLORORDER;COLORORDER];


%this is sort of set by hand to get everything fitting
% by changing the size of the axes in the plot
% works for the 3 line title 
%SQUEEZE = [1 1 0.93 0.88];
SQUEEZE = [1 1 1 1];

% children, like title, xlabel, etc of axes
% this will be applied if this function called after plot made
CFONTSIZE = 10;

HA = get(HF,'children');
for jj = 1:length(HA)
    set(HA(jj),'box','on');
    set(HA(jj),'fontsize',FONTSIZE,'fontname',FONTNAME);
    AA = get(HA(jj),'position');
    set(HA(jj),'position',AA.*SQUEEZE);

    % COMMENT this section out (or at least the line stuff)
    % if you do not want the line sizes etc to be changed

    HL = get(HA(jj),'children');
    HT = get(HA(jj),'title');
    HX = get(HA(jj),'xlabel');
    HY = get(HA(jj),'ylabel');

    HL = [HL;HT;HX;HY];

	    for ii = 1:length(HL);

		    TYPE = get(HL(ii),'type');

		    if strcmp(TYPE,'line');
			    set(HL(ii),'linewidth',LINEWIDTH);
			    set(HL(ii),'color',COLORORDER(ii));

		    elseif strcmp(TYPE,'text');
			    set(HL(ii),'fontname',FONTNAME,'fontsize',CFONTSIZE);

		    end

	    end

    end