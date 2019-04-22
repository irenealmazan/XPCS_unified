%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if any(WhatToDo==4) %% make summed images 
%	ImageJ = 1;  % this must be set and be consistent earlier! 
%	ImageJ indexing (start (0,0)) in ROI indexing
	 [SIZE1,SIZE2,SIZE3]=size(IInormb);  % in case whole scan not complete
	 [SoY,SoX] = slicesumrois(IInormb,ROIS,ImageJ);%(1,:));
%	Note - This does not 'normalize' the slices to the number of pixels summed
%			However, after using function, use SoY.image{i}./SoY.norm{i} 
%				for ROI {i} 

% N_ROIS_SUM vector with the ROIS to use, SoX(Y)FLAG whether to sum in X or Y for that ROI
if length(SoXFLAG)~=length(N_ROIS_SUM);
	SoXFLAG = ones(1,length(N_ROIS_SUM));
end
if length(SoYFLAG)~=length(N_ROIS_SUM);
	SoYFLAG = ones(1,length(N_ROIS_SUM));
end	

iii = 0;
for ii=N_ROIS_SUM; 
  iii=iii+1;
		
  
  if SoXFLAG(iii)
	figure;clf
	set(gcf,'Name',['SoX ROIs #' int2str(ii)]);

	 if ~LOGFLAG;
		HS = imagesc(Xsteps,SoX.ndx{ii},SoX.images{ii});
%		HS = surf(Xsteps,SoX.ndx{ii},SoX.images{ii}, SoX.images{ii}));
	 else
		HS = imagesc(Xsteps,SoX.ndx{ii},log10(SoX.images{ii}));
%		HS = surf(Xsteps,SoX.ndx{ii},log10(SoX.images{ii}),log10(SoY.images{ii}));
	 end

	 HPC.SoX_ROIS(ii) = HS;
		shading flat; view(0,90); colormap(COLORMAP)
		xlabel(SCNXLABEL);ylabel(YROW);
		title(char(TITLEstr1,[DOCUInt, ' summed over XCOL in ROI #' int2str(ii)]));
	 set(gca,'Xlim',[min(Xsteps) max(Xsteps)]);
	 
		% Now YLim is simply Y on the graph
		%	if ~isempty(AXISdet); set(gca,'Ylim',AXISdet(3:4));end
		if diff(ROIS(ii,[3:4]))   % if non zero do this
			set(gca,'Ylim', ROIS(ii,[3:4])); 
	 	else % accommodate line ROI that user summs in non-useful dirction
			set(gca,'Ylim', ROIS(ii,[3])+[-5 5]);
		end

	 			
	plot_adjust(gca,[0 1]);
	 
		
		if ~isempty(lastframes)
			HYline = makeyline(Xsteps(lastframes+ImageJ),'b',gca);  % and makes a straight line
		end

  end
	
  if SoYFLAG(iii)  % don't plot if SoY flag 0
	figure;clf
	set(gcf,'Name',['SoY ROIs #' int2str(ii)]');
	
	 if ~LOGFLAG;  % surf(X,Y,Z,Color)
%	 HS = surf(Xsteps,SoY.ndx{ii},((SoY.images{ii})),((SoY.images{ii})));	
	 HS = imagesc(Xsteps,SoY.ndx{ii},((SoY.images{ii})));		
	 else
%	 HS = surf(Xsteps,SoY.ndx{ii},log10(fixlogimage(SoY.images{ii})),log10(fixlogimage(SoY.images{ii})));
	 HS = imagesc(Xsteps,SoY.ndx{ii},log10(SoY.images{ii}));
	 end

	 HPC.SoY_ROIS(ii) = HS;
	 shading flat; view(0,90); colormap(COLORMAP)
	 xlabel(SCNXLABEL);ylabel(XCOL);
	 title(char(TITLEstr1,[DOCUInt ' summed over YROW in ROI #' int2str(ii)]));
	 set(gca,'Xlim',[min(Xsteps) max(Xsteps)]);

		% Now YLim is simply Y on the graph
		%	if ~isempty(AXISdet); set(gca,'Ylim',AXISdet(3:4));end
		if diff(ROIS(ii,[1:2]))   % if non zero do this
			set(gca,'Ylim', ROIS(ii,[1:2])); 
	 	else % accommodate line ROI that user summs in non-useful dirction
			set(gca,'Ylim', ROIS(ii,[1])+[-5 5]);
		end



	plot_adjust(gca,[0 1]);

		if ~isempty(lastframes)
			HYline = makeyline(Xsteps(lastframes+1),'b',gca);  % and makes a straight line
		end

  end
end

	if  1    %   plots that show line plots as function of nu or del
		% Note - if any SoX (SoY) requested  on any ROIs, the all ROIS will done
		if any(SoXFLAG)
		figure; clf
		set(gcf,'Name','SoX and Sum Points ');
		HL = semilogy(SoX.ndx{1},sum(SoX.images{1}'));
					% in event only one point because of dimensionality
					if size(sum(SoX.images{1}'))<2;HL(1).Marker='o';end
		
			for jj = 2:length(SoX.images);
				HL(jj) = line(SoX.ndx{jj},sum(SoX.images{jj}'));
					% in event only one point because of dimensionality
					if size(sum(SoX.images{jj}'))<2;HL(jj).Marker='o';end
			end

		xlabel(YROW),ylabel('Int (arb)'),
		title(char(TITLEstr1,['summed over XCOLS in all ROI and ',int2str(SIZE3),' scan pts']));
		plot_adjust(gca,0); % 
			for ii=1:length(ROIS(:,1));
				set(HL(ii),'Color',COLORORDER(ii,:),'LineWidth',2,'DisplayName',['ROI #',int2str(ii)]);
			end
				Hlegend = legend('show','AutoUpdate','off');

		
		
		HLine.SoX_SPts = HL;clear HL;
	
		end
	
		if any(SoYFLAG)	
		figure; clf
		set(gcf,'Name','SoY and Sum Points ');
		HL = semilogy(SoY.ndx{1},sum(SoY.images{1}'));
							% in event only one point because of dimensionality
					if size(sum(SoY.images{1}'))<2;HL(1).Marker='o';end
			for jj = 2:length(SoY.images);
				HL(jj) = line(SoY.ndx{jj},sum(SoY.images{jj}'));
					% in event only one point because of dimensionality
					if size(sum(SoY.images{jj}'))<2;HL(jj).Marker='o';end
			end

		xlabel(XCOL),ylabel('Int (arb)'),
		title(char(TITLEstr1,['summed over YROWS in all ROI and ',int2str(SIZE3),' scan pts']));
		plot_adjust(gca,0); %
			for ii=1:length(ROIS(:,1));
				set(HL(ii),'Color',COLORORDER(ii,:),'LineWidth',2,'DisplayName',['ROI #',int2str(ii)]);
			end
				Hlegend = legend('show','AutoUpdate','off');
	 
		HLine.SoY_SPts = HL;clear HL;
		end


		if ~isempty(POINTSUMS)   % sum over points sets in POINTSUMS, one plot per ROI
	
		jjj=0; 
	
			for jj=N_ROIS_SUM;
				jjj=jjj+1;
				
				idROI = int2str(jj);
				if SoXFLAG(jjj)
		
				figure; clf;
				set(gcf,'Name',['SoX ROI # ' int2str(jj) ' and Sum select points']);
					for ii=1:length(POINTSUMS(:,1));
					Ni = [POINTSUMS(ii,1): POINTSUMS(ii,2)]+ImageJ;  %use as matlab
%					HL(ii) = line(SoX.ndx{1},sum(SoX.images{1}(:,Ni)'));
					HL(ii) = line(SoX.ndx{jj},sum(SoX.images{jj}(:,Ni)'./length(Ni)));
					end
		
		
				set(gca,'Yscale','log')
				xlabel(YROW),ylabel('Int (arb)'),
				title(char(TITLEstr1,['# ' int2str(jj) ' ROI summed over X, then summed over Spec Point Range']));
				plot_adjust(gca,0); %
					legend(addnames2matrix('sumX between [', int2str(POINTSUMS),']/Npts'))
				HLine.SoX_SoSelectPt = HL;clear HL;
				set(HLine.SoX_SoSelectPt,'linewidth',1.5);
				end
		
				if SoYFLAG(jjj)

				figure; clf;
				set(gcf,'Name',['SoY ROI # ' int2str(jj) ' and Sum selected points']);
					for ii=1:length(POINTSUMS(:,1));
						Ni = [POINTSUMS(ii,1): POINTSUMS(ii,2)]+ImageJ;
%						HL(ii) = line(SoY.ndx{1},sum(SoY.images{1}(:,Ni)'));
						HL(ii) = line(SoY.ndx{jj},sum(SoY.images{jj}(:,Ni)'./length(Ni)));
					end
				set(gca,'Yscale','log')
				xlabel(XCOL),ylabel('Int (arb)'),
				title(char(TITLEstr1,['# ' int2str(jj) ' ROI summed over Y, then summed over Spec Point Range']));
				plot_adjust(gca,0); %
					legend(addnames2matrix('sumY between [', int2str(POINTSUMS),']/Npts'));
				HLine.SoY_SoSelectPt = HL;clear HL;
				set(HLine.SoY_SoSelectPt,'linewidth',1.5);
				end
			end   % end of if jj % jj=N_ROIS_SUM;
	
		end  % sum over points sets in POINTSUMS
	end %%   plots that show line plots as function of nu or del
% end
