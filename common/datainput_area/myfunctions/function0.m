%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if any(WhatToDo==0); %% plot a 'single' image from the scan

	if ~LOGFLAG
		figure;clf
		set(gcf,'Name',['Lin Image',THRESHusedoc,' Spec # ' int2str(SINGLE)]);
		HS = imagesc(XCOLpts,YROWpts,(IInormb(:,:,SINGLE+ImageJ)));
		shading interp
		HSurf.IInormb = HS;
		
			
			if ~ASPECTflag
			axis equal;
			else
			HA = gca;
			daspect(gca,ASPECT);   %this might not be in earlier version
			HA = gca;,
			set(HA,'DataAspectRatioMode','manual','DataAspectRatio',[ASPECT]);
			end
			
			if ~isempty(AXISdet); axis(AXISdet);end
			if ~isempty(CLIMlin); 
				set(gca,'clim',CLIMlin);
				DOCUclim = mDOCUclim(CLIMlin);
			end;
	else
		figure;clf
		set(gcf,'Name',['Log10 Image',THRESHusedoc,' Spec #' int2str(SINGLE)]);  % assumes SINGLE in spec/ImageJ
		HS = imagesc(XCOLpts,YROWpts,log10((IInormb(:,:,SINGLE+ImageJ))));
		shading flat
		HSurf.IInormb = HS;

			
			if ~ASPECTflag
			axis equal;
			else
			HA = gca;
			daspect(gca,ASPECT);   %this might not be in earlier version
			HA = gca;,
			set(HA,'DataAspectRatioMode','manual','DataAspectRatio',[ASPECT]);
			end
			
			
			if ~isempty(AXISdet), axis(AXISdet);end
			if ~isempty(CLIMlog); 
				set(gca,'clim',CLIMlog);
				DOCUclim = mDOCUclim(CLIMlog);
			end;
	end
	
	DOCUlast = ['(Spec pt #',int2str(SINGLE),') ', DOCUInt, DOCUclim];
	title(char(TITLEstr1,DOCUlast,pfilename(INFOstr)));
	xlabel(XCOL);ylabel(YROW);
	%prettyplot(gca,[2 2 4 3.8]);
	plot_adjust(gca,[0 1],[],[2 2 4 3.8]);
	colormap(COLORMAP);
	
% end	
