%%%%%%%%%%%%%%%%  HELPER FUNCTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DOCUclim = mDOCUclim(HA);
% to have a method to update clim for legend/titles

	if nargin<1;HA=gca;end
	if isempty(HA)
		DOCUclim = [];
	elseif ~(length(HA)==2);
		CLIM=get(gca,'clim');
		DOCUclim = [' : clim[', num2str(CLIM), ']'];
	else length(HA)==2;
		CLIM=HA;
		DOCUclim = [' : clim[', num2str(CLIM), ']'];
	end
end
	


