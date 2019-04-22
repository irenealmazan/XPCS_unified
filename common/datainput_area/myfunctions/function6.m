%if any(WhatToDo==6);   %to make it automatic for summs   % ? is Xsteps missing one point
	if isempty(POINTSUMS); 
		POINTS=[1 : length(Xsteps)]-ImageJ;
	else 
		POINTS = POINTSUMS;
	end

for ii = 1:length(POINTS(:,1))
	Ni = [POINTS(ii,1): POINTS(ii,end)]+ImageJ;
	NP = length(Ni);
	POINTstr = ['[ ' int2str(Ni(1)-ImageJ) ':' int2str(Ni(end)-ImageJ) ']'];
	
	figure;clf;
	set(gcf,'Name','Image summed over points/Num Points');	
	if ~LOGFLAG

		HS = imagesc(XCOLpts,YROWpts,sum(Iwinnorm(:,:,Ni),3)./(NP));
		colormap(COLORMAP)
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
			else 
				CLIMtemp = get(gca,'clim');
				DOCUclim = mDOCUclim(CLIMtemp);
			end;
			DOCUInt6 = [DOCUInt, ' SUMMED over SPEC scan pts ',POINTstr];
	else
	
		HS = imagesc(XCOLpts,YROWpts,log10(sum(IInormb(:,:,Ni),3)./(NP)));
		colormap(COLORMAP)	
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
			else 
				CLIMtemp = get(gca,'clim');
				DOCUclim = mDOCUclim(CLIMtemp);
			end;
			DOCUInt6 = [DOCUInt, ' SUMMED over SPEC scan pts ', POINTstr];
	end
	
	shading flat;
	
	if FLATFIELDFLAG
		DOCUInt = [DOCUInt,' [FF]'];
	end
	
	DOCUlast = char(DOCUInt6,pfilename(INFOstr));
	title(char(TITLEstr1,DOCUlast));xlabel(XCOL);ylabel(YROW);
	plot_adjust(gca,0,[],[2 2 4 3.8]); %
		% only put ROIs on if we did something like option 2
		if any(WhatToDo==2);
			if exist('COLORORDER','var')%length(ROIS(:,1))<length(COLORORDER(:,1))
				showrois(ROIS,gcf,1.5,COLORORDER);
			else
				[HROI,COLORORDER]=showrois(ROIS,gcf,1.5);
			end
		end
	
end

%end
