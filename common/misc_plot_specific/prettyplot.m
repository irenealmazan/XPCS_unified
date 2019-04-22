function [PRESET] = prettyplot(HA,PAPERPOSITION,SQUEEZEflag);
% [PRESET] = prettyplot(HA,PAPERPOSITION,SQUEEZEflag) Set the papersize, and the fontsizes and types 
% 
%	OUTPUT
%		[PRESET] is mostly just a flag
%	INPUTs are optional as these are mostly assumed to be hardwired for particular set of figures
%		HA 				(Axed graphic handle) Default assumes HA=gca
%		PAPERPOSITION	(4-vector)   position (X,Y) in inches from lower left corner, then width 
%						Default [1 1 3.75 3.5]   i.e., [Xpos,Ypos,Xwidth,Ywidth]
%		SQUEEZEflag		(scalar) =1 (apply squeeze, Default) = 0 (do not apply squeeze)
%						As we typically have multiline title this helps squeezes the axes
%						Note - it will keep squeezing each time applied
%
%						check interior of program
%						SQUEEZE (4 element vector) often gets tweaked 
%						SQUEEZE=[1 1 1 1] will have no effect, 
%						SQUEEZE=[1 1 0.95 .9] has typically been applied to help the multiple line title not get shoved off
%
%	  Use plot_adjust instead 
%			plot_adjust(HA,0) to get the same effect as prettyplot

if nargin<3;SQUEEZEflag=1;end
if nargin<2;PAPERPOSITION = [2 2 3.75 3.5];end
if nargin<1;HA = gca;end

colorflag = [0 1];  % do not alter colors of lines, [0] but do alter linewidths [1]
SQUEEZE = [1 1 1 1];
% The following have at various times been used so multiline title fits onto prints
% SQUEEZE = [1 1 0.93 .88];
% SQUEEZE = [1 1 0.95 0.9];

% basically, plot_adjust and do not want to keep up two programs
% prettyplot was basically one that fussed with the fonts and print size,
%  and did not change colors of lines 
% although I usually expected that it might make them a bit wider from defaults
[HL,HText] = plot_adjust(HA,colorflag,SQUEEZE,PAPERPOSITION);

