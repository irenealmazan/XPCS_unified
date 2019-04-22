function [ROIsconv] = ROISCEN2ROIS(ROIS,CEN)
% [ROIsconv] = ROIS2CENROIS(ROIS,[newXCEN,newYCEN)  shift ROIS from one center to another
%  input ROIS, CENS = [newXCEN newYCEN;oldXCEN oldYCEN]
if nargin<2;help ROISCEN2ROIS;return;end % 


%
	ROIsconv = ...
		[([ROIS(:,1) ROIS(:,2)] - diff(CEN(:,1))) ([ROIS(:,3) ROIS(:,4)]-diff(CEN(:,2)))];
