%%% note if the ROI has NaN in any of the pixels, this will mess up and give NaN
%%% need to fix sumrois to fix that
% if any(WhatToDo==5)  %% make summed lines of ROI's on images
	figure;clf
	set(gcf,'Name','ROIs as counters');
	 
	 [sumROIS,sumROISnormed] = sumrois(IInormb,ROIS,ImageJ);
	 %  sumROIS sums over the ROI, sumROISnormed divides by the number of pixels
%	 [sumROIS,sumROISnormed] = sumrois(IInormb,ROIS,ImageJ);
	 %	Isum	Array  : each column for each ROI. Length of vectors (rows) num of images
	 %HL = semilogy(Xsteps,sumROIS);
	 HL = plot(Xsteps,sumROISnormed); 
		YLABEL = '[Inorm(summed over ROI)/ROIsize]';
%	 HL = plot(Xsteps,log10(sumROISnormed)); YLABEL = 'log10[Inorm(summed over ROI)/ROIsize]';


	 xlabel(SCNXLABEL);ylabel(YLABEL);
	 title(char(TITLEstr1,pfilename(INFOstr),YLABEL));
	 plot_adjust(gca,0); %
% mocvd
%	 makeyline(Xsteps(lastframes+ImageJ),'b',gca);   % want colored lines on this plot, not white	 
	 for ii=1:length(ROIS(:,1));
		set(HL(ii),'Color',COLORORDER(ii,:),'LineWidth',2,'DisplayName',['ROI #',int2str(ii)]);
	 end

	 Hlegend = legend({HL.DisplayName},'AutoUpdate','off');   % why does it work here and not in 4?
	 makeyline(Xsteps(lastframes+ImageJ),'b',gca);   % want colored lines on this plot, not white	 
	 
	 
% end
