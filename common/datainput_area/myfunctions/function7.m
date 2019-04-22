% if any(WhatToDo==7);  % plot flatfield FLATFIELDFLAG

	figure
	set(gcf,'Name',['Flat Field']);
	HS = imagesc(XCOLpts,YROWpts,FLATFIELD);
	shading flat
	axis equal;
	colormap(COLORMAP)
    if ~isempty(AXISdet); axis(AXISdet);end
	title(['Flat Field: from ' pfilename(FlatField.filenameroot) ' [' FLATFIELDinfo ']']);
	xlabel(XCOL);ylabel(YROW);
	plot_adjust(gca,0,[],[2 2 4 3.8]); %
    %set(gca,'position',[.17 .1 .65 .68]);

			if exist('COLORORDER','var')%length(ROIS(:,1))<length(COLORORDER(:,1))
				showrois(ROIS,gcf,1,COLORORDER);
			else
				[HROI,COLORORDER]=showrois(ROIS,gcf,1);
			end	


% end