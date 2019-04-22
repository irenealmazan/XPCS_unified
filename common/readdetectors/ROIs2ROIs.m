function [ROIsconv] = ROIS2ROIS(ROIS,direction)
% [ROIsconv] = ROIS2ROIS(ROIS,direction)  change the ROIS between [STARTpix DELTApix] to [STARTpix ENDpix]
%	(n.b. this isn't changing whether indexing scheme starts at 1 or 0!)
%
% Output 
%		ROIsconv (matrix) of converted ones 
% Input
%		ROIs (matrix) of initial values  (assumed X X Y Y) each row new ROI
%		direction = 1 or 0; 1 is default
%				1  convert [STARTpix DELTApix] to [STARTpix ENDpix]  (to work in my programs)
%				0  convert [STARTpix ENDpix] to [STARTpix DELTApix]
%
% The matlab functions of CT happen to be set up as STARTpix ENDpix
% But the ROI's as put in the EPICS medm (and from spec file) are 
%		Aversion				Bversion
%  STARTpix DELTApix <-->  STARTpix STARTpix+DELTApix-1
%
%
% Example - in [STARTpix DELTApix] an ROI might look like
%	[Xstart Xwidth Ystart Ywidth]
%	[100 5 200 6]	5 X by 6 Y pixel ROI with lower corners at 100,200 [X,Y]
%
% Is equivalent to [STARTpix ENDpix] ROI that looks like
%	[Xstart Xend Ystart Yend]
%	[100 104 200 205]	5 X by 6 Y pixel ROI with lower corners at 100,200 [X,Y]
%

%
%  Note - a way to remember the extra -1 or +1 is to think about
%		if ROI were single pixel wide, e.g., (or at least I think that we'd put a
%			1 in the epics MEDM DELTApix if the ROI were 1 pixel wide
%  5 1  <-->  5 5
%  e.g., (at least I think that we'd put a 1 in the epics MEDM DELTApix 
%           if we wanted an ROI that were 1 pixel wide

if nargin<2; direction=1;end % default [STARTpix DELTApix] --> [STARTpix ENDpix] 

if direction~=0;
	ROIsconv = ...
		[ROIS(:,1) ROIS(:,2)+ROIS(:,1)-1 ROIS(:,3) ROIS(:,4)+ROIS(:,3)-1];
else 

	ROIsconv = ...
		[ROIS(:,1) ROIS(:,2)-ROIS(:,1)+1 ROIS(:,3) ROIS(:,4)-ROIS(:,3)+1];
end
