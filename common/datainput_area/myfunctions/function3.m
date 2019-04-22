% if any(WhatToDo==3)  %% make movie  (frame must be 560x413?)
	figure;clf;figure(gcf);
	set(gcf,'Name',['MOVIE ' THRESHusedoc]);
	% set up initial size and image and climits  (Default SINGLE is 0 (spec point + ImageJ = matlab point)
	%	
	% Note that video get upset if size changes of image some of this is to 
	% try to make sure this doesn't happen but it still happens and I 
	% don't know how to fix at times except to keep trying again (magic)
	clear HS;
	if ~LOGFLAG
		HS = imagesc(XCOLpts,YROWpts,((IInormb(:,:,SINGLE+ImageJ))));  
	else
		HS = imagesc(XCOLpts,YROWpts,log10(IInormb(:,:,SINGLE+ImageJ)));  
	end
	
	colormap(COLORMAP)
	shading flat
		HT = title(...
			char([TITLEstrshort],...
			['Spec Pt # [ 00000 ] : step =  [0000000]',mDOCUclim(CLIMlog)]));	

		if ~ASPECTflag
		axis equal;
		else
			HA = gca;
			daspect(gca,ASPECT);   %this might not be in earlier version
			HA = gca;,
			set(HA,'DataAspectRatioMode','manual','DataAspectRatio',[ASPECT]);
		end
			
		
		xlabel(XCOL);ylabel(YROW);
		
		if MOVIEROIFLAG
			MAA = ROIS(1,:);
		else
			if ~isempty(AXISdet);
				MAA = AXISdet;
			else
				MAA = [XCOLpts(1) XCOLpts(end) YROWpts(1) YROWpts(end)];
			end
		end	

		set(gca,'xlim',MAA([1 2]),'ylim',MAA([3 4]));
		
		if ~isempty(CLIMlog) & LOGFLAG; 
			set(gca,'clim',CLIMlog);
		elseif ~isempty(CLIM) & ~LOGFLAG;
			set(gca,'clim',CLIMF);
		end;
		plot_adjust(gca,[0 1]);
		
		if any(WhatToDo==2);
			if exist('COLORORDER','var')%length(ROIS(:,1))<length(COLORORDER(:,1))
				showrois(ROIS,gcf,0.5,COLORORDER);
			else
				[HROI,COLORORDER]=showrois(ROIS,gcf,0.5);
			end
		end		
		

	% Prepare the new file.
		if exist('M1.avi')==2; delete('M1.avi'); end  
		vidObj = VideoWriter('M1.avi');
		vidObj.FrameRate = 5;		% default 30 fps


		open(vidObj);


%	for ii=1:length(timestampX)  % CT
	for ii=1:length(Xsteps);
		PTstr = int2str(ii-1+10000);PTstr = PTstr(2:end);
		
		if ~LOGFLAG
			set(HS,'CData',((IInormb(:,:,ii))));
			if ~isempty(CLIMlin); set(gca,'clim',CLIMlin);end;
		else
			set(HS,'CData',(log10(IInormb(:,:,ii))));
			if ~isempty(CLIMlog); set(gca,'clim',CLIMlog);end;
		end		
		set(HT,'String',...
			char([TITLEstrshort],...
			['Spec Pt # [',PTstr,'] : step = [',num2str(Xsteps(ii)),']',mDOCUclim(gca)]));
		currFrame = getframe(gcf);
		writeVideo(vidObj,currFrame);
	end
	
	% have short pause at end of movie while sitting on the end point
	for ii=1:3;
		currFrame = getframe(gcf);
		writeVideo(vidObj,currFrame);
	end
	
	disp('*');
	%Mov_n160208a_intval=M;
	close(vidObj);
	
% end
