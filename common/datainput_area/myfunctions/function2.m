
%if any(WhatToDo==2)   % show the ROIS lines on clean figure that looks like the previous figure via AXISdet or not

	[HROI,COLORORDER] = showrois(ROIS,gcf);  % first put the ROIS on the current figure
		
		figure;clf;  % then makes a new clean figure to show them
		set(gcf,'Name','Enumeration of ROIs colors');

%		prettyplot(gca);    % put here or it clobbers the showROIS colororder 		
		showrois(ROIS,gcf,2,COLORORDER);
		title(char('Enumeration of the ROIs and their colors used for',...
			[pfilename(sdata.filename),' #', int2str(sdata.SCNNUM),' : ', sdata.SCNDATE]));

		axis equal;
		if ~isempty(AXISdet); 
			axis(AXISdet);
		else
			axis([min(XCOLpts) max(XCOLpts) min(YROWpts) max(YROWpts)]);
		end

			LEGEND = [char([ones(length(HROI),1)*'ROI #']) int2str([1:length(HROI)]')];
			legend(LEGEND);
			xlabel(XCOL);ylabel(YROW);  % must have done the WhattoDo =0
			
			plot_adjust(gca,[0 1]);
		
%end
