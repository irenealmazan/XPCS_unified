classdef Function_display_sos
    % This library contains all the functions to analyze the simulations
    properties(Constant)
    end
    
    
    methods(Static)
        
        function [hl] = show_ranges_averaged(QX,QY,iyaoff,iyahw,ixaoff,ixahw)
            
            if numel(QX)~= numel(QY)
                QX = QY; 
           
            end
                
            
            % Show ranges averaged
            for kk = 1:numel(iyaoff)
                hl = line(QX,(iyaoff(kk)+iyahw)*ones(size(QY)));
                set(hl,'Linestyle','--','Color','w','LineWidth',4.0);
                hl = line(QX,(iyaoff(kk)-iyahw)*ones(size(QY)));
                set(hl,'Linestyle','--','Color','w','LineWidth',4.0);
                
            end
            
            for kk = 1:numel(ixaoff)
               hl = line((ixaoff(kk)+ixahw)*ones(size(QX)),QY);
                set(hl,'Linestyle','--','Color','w','LineWidth',4.0);
                hl = line((ixaoff(kk)-ixahw)*ones(size(QX)),QY);
                set(hl,'Linestyle','--','Color','w','LineWidth',4.0);
            end
        end
        
        function [legend_str] = show_multiple_offsets(QX,QY,iyaoff,iycen,iyahw,ixaoff,ixcen,ixahw,tauML,color)
            
            if numel(iyaoff) > 1
                for kk = 1:numel(iyaoff)
                    hl = line((QX-ixaoff),mean(tauML(iycen+iyaoff(kk)+[-iyahw:iyahw],:),1),'Color',color(kk));
                    set(hl,'LineStyle','none','Marker','o');
                    set(hl,'LineStyle','none','Marker','o');
                    legend_str{kk} = ['offset in Q_y = ' num2str(iyaoff(kk))];
                end
            elseif numel(ixaoff)>1
                for kk = 1:numel(iyaoff)
                    hl = line(QY,mean(tauML(:,ixcen+ixaoff(kk)+[-ixahw:ixahw]),2),'Color',color(kk));
                    set(hl,'LineStyle','none','Marker','o');
                    legend_str{kk} = ['offset in Q_x = ' num2str(ixaoff(kk))];
                end
            end
            
        end
        
        function [] = make_figure_2D(title_figure,map2D,field_to_plot,struct_offsets,iCTR,figNum,varargin)
            
             p = inputParser;
            
            addRequired(p,'title_figure');
            addRequired(p,'map2D');
            addRequired(p,'field_to_plot');
            addRequired(p,'struct_offsets',@isstruct);
            addRequired(p,'iCTR',@isnumeric);
            addRequired(p,'figNum',@isnumeric);
            
            addParameter(p,'POSITION',[100   100   400   300],@isnumeric);
            addParameter(p,'PAPERPOSITION',[1     1     4     3],@isnumeric);
            addParameter(p,'PAPERSIZE',[ 8.5000   11.0000],@isnumeric);           
            addParameter(p,'CLIM',[6.1971    9.3297],@isnumeric);
            addParameter(p,'plot_averages',1,@isnumeric);
            addParameter(p,'ZSCALE','linear',@ischar);
            
            parse(p,title_figure,map2D,field_to_plot,struct_offsets,iCTR,figNum,varargin{:});
            
            Nsubplot = numel(map2D);
            
            figure(figNum);
            clf;
            
            for kk = 1:Nsubplot
                
                subplot(1,Nsubplot,kk);
                hl = pcolor(struct_offsets.QX,struct_offsets.QY,log10(map2D(kk).(field_to_plot)));
                set(gca,'Clim',p.Results.CLIM);

                set(gca,'ZScale',p.Results.ZSCALE);
                %axes('Box','on');
                shading flat;
                axis image;
                xlabel('Q_X (pixels)');
                ylabel('Q_Y (pixels)');
                colorbar;
                
                % Show ranges averaged                
                if p.Results.plot_averages
                    Function_display_sos.show_ranges_averaged(struct_offsets.QX,struct_offsets.QY,struct_offsets.yaoff_array,struct_offsets.iyahw,struct_offsets.xaoff_array,struct_offsets.ixahw);
                end
                
                % show CTR
                hl = line(struct_offsets.nsteps*iCTR,zeros(size(iCTR)));
                set(hl,'Linestyle','none','Marker','.','Color','r');
                
                title([title_figure(kk,:)]);
                
            end
            
            set(gcf,'Position',p.Results.POSITION);
            set(gcf,'PaperPosition',p.Results.PAPERPOSITION);
            set(gcf,'PaperSize',p.Results.PAPERSIZE);
            %title([runname_title ' Correlation Time (ML) in log scale']);
            
        end
        
        function [] = make_figures_tauML1D(title_figure,Q,tauML_1D_struct,field_to_plot,figNum,varargin)
            
             p = inputParser;
            
            addRequired(p,'title_figure');
            addRequired(p,'Q');
            addRequired(p,'tauML_1D_struct');
            addRequired(p,'field_to_plot');
            addRequired(p,'figNum',@isnumeric);
            
            addParameter(p,'tauMLth',[],@isnumeric);
            addParameter(p,'POSITION',[100   100   400   300],@isnumeric);
            addParameter(p,'iCTR',[],@isnumeric);
            addParameter(p,'ixaoff',0,@isnumeric);
            addParameter(p,'iyaoff',0,@isnumeric);
            addParameter(p,'nsteps',@isnumeric);
            addParameter(p,'PAPERPOSITION',[1     1     4     3],@isnumeric);            
            addParameter(p,'YLIM',[1e6    1e10],@isnumeric);
            addParameter(p,'XLABEL','',@ischar);
            addParameter(p,'XSCALE','linear',@ischar);
            addParameter(p,'YSCALE','linear',@ischar);
            addParameter(p,'FONTSIZE',15,@isnumeric);

            parse(p,title_figure,Q,tauML_1D_struct,field_to_plot,figNum,varargin{:});
            

            % Display tauML integrated over rows:
            figure(figNum);
            clf;
            
            axes('Box','on');

            hold on;
           
            for kk = 1:numel(tauML_1D_struct)
                % hl = line(struct_offsets.QX,mean(tauML(range_rows,:),1));
                hl = line(Q,tauML_1D_struct(kk).(field_to_plot));
                
                if isfield(tauML_1D_struct(kk),'Color')
                    set(hl,'LineStyle','none','Marker','o','Color',tauML_1D_struct(kk).Color);                    
                else
                    set(hl,'LineStyle','none','Marker','o','Color',[0.0 0.45 0.74]);                   
                end
            end
            set(gca,'Ylim',p.Results.YLIM);
            pa = axis;
            
            if ~isempty(p.Results.tauMLth)
                hl = line(Q,p.Results.tauMLth);
                set(hl,'LineStyle','-','Color','m','LineWidth',2.0);
                axis(pa);
            end
            
            % Show CTR positions
            if ~isempty(p.Results.iCTR)
                for ii = p.Results.iCTR
                    hl = line(p.Results.nsteps*ii*[1 1]-p.Results.ixaoff,pa(3:4));
                    set(hl,'LineStyle','--','Color','r');
                end
            end
            
            set(gca,'YScale',p.Results.YSCALE,'XScale',p.Results.XSCALE);
           
            xlabel(p.Results.XLABEL);
            ylabel('Mean Correlation Time (s)');
            title([title_figure]);
            %legend(['Offset Q_y = ' num2str(struct_offsets.iyaoff)]);
            set(gcf,'Position',p.Results.POSITION);
            set(gcf,'PaperPosition',p.Results.PAPERPOSITION);
            set(gca,'FontSize',p.Results.FONTSIZE);
            
            
        end
        
        function [] = make_figures_tauML_1D_log(title_figure,tauML,struct_offsets,iCTR,figNum,varargin)
            
             p = inputParser;
            
            addRequired(p,'title_figure');
            addRequired(p,'tauML');
            addRequired(p,'struct_offsets',@isstruct);
            addRequired(p,'iCTR',@isnumeric);
            addRequired(p,'figNum',@isnumeric);
            
            addParameter(p,'tauML_th_QX',[],@isnumeric);
            addParameter(p,'tauML_th_QY',[],@isnumeric);
            addParameter(p,'POSITION',[100   100   400   300],@isnumeric);
            addParameter(p,'PAPERPOSITION',[1     1     4     3],@isnumeric);            
            addParameter(p,'YLIM',[1e6    1e10],@isnumeric);
            
            parse(p,title_figure,tauML,struct_offsets,iCTR,figNum,varargin{:});
            
            
            
            range_rows = struct_offsets.iycen+struct_offsets.iyaoff+[-struct_offsets.iyahw:struct_offsets.iyahw];
            tauML_Qx = mean(tauML(range_rows,:),1);
            
            Qx_positive = struct_offsets.QX(struct_offsets.QX>0);
            tauML_Qx_positive = tauML_Qx(struct_offsets.QX>0);
            tauML_th_QX_positive = tauML_th_QX(struct_offsets.QX>0);
            
            Qx_negative = struct_offsets.QX(struct_offsets.QX<=0);
            tauML_Qx_negative = tauML_Qx(struct_offsets.QX<=0);
            tauML_th_QX_negative = tauML_th_QX(struct_offsets.QX<=0);
            
            
            range_cols = struct_offsets.ixcen+struct_offsets.ixaoff+[-struct_offsets.ixahw:struct_offsets.ixahw];
            tauML_Qy = mean(tauML(:,range_cols),2);
            
            Qy_positive = struct_offsets.QY(struct_offsets.QY>0);
            tauML_Qy_positive = tauML_Qy(struct_offsets.QY>0);
            tauML_th_QY_positive = tauML_th_QY(struct_offsets.QY>0);
            
            Qy_negative = struct_offsets.QY(struct_offsets.QY<=0);
            tauML_Qy_negative = tauML_Qy(struct_offsets.QY<=0);
            tauML_th_QY_negative = tauML_th_QY(struct_offsets.QY<=0);
            
            % Display tauML integrated over rows:
            figure(figNum);
            clf;
            
            axes('Box','on');

            hold on;
            
            subplot(121)
            Function_display_sos.make_subplot_log(Qx_positive , tauML_Qx_positive,...
                'tauML_th',tauML_th_QX_positive,'iCTR',iCTR,...
                'nsteps',struct_offsets.nsteps,'XLABEL','Q_x (pixels)',...
                'TITLE','positive Qx');
            
            subplot(122)
            Function_display_sos.make_subplot_log(Qx_negative , tauML_Qx_negative,...
                'tauML_th',tauML_th_QX_negative,'iCTR',iCTR,...
                'nsteps',struct_offsets.nsteps,'XLABEL','Q_x (pixels)',...
                'TITLE','negative Qx');
            
            figure(figNum+2);
            clf;
            
            axes('Box','on');

            hold on;
            
            subplot(121)
            Function_display_sos.make_subplot_log(Qy_positive , tauML_Qy_positive,...
                'tauML_th',tauML_th_QY_positive,...
                'XLABEL','Q_y (pixels)',...
                'TITLE','positive Qy');
            
            subplot(122)
             Function_display_sos.make_subplot_log(Qy_negative , tauML_Qy_negative,...
                'tauML_th',tauML_th_QY_negative,...
                'XLABEL','Q_y (pixels)',...
                'TITLE','negative Qy');
            
            
            
            %{
            hl = line(Qx_positive , tauML_Qx_positive);
            
            set(hl,'LineStyle','none','Marker','o','Color',[0.0 0.45 0.74]);
            pa = axis;
            
            if ~isempty(p.Results.tauML_th_QX)
                hl = line(Qx_positive,tauML_th_QX_positive);
                set(hl,'LineStyle','-','Color','m','LineWidth',2.0);
                axis(pa);
            end
            
            % Show CTR positions
            for ii = iCTR
                hl = line(struct_offsets.nsteps*ii*[1 1],pa(3:4));
                set(hl,'LineStyle','--','Color','r');
            end
                        
            set(gca,'Ylim',p.Results.YLIM);
            xlabel('Q_X (pixels)');
            ylabel('Mean Correlation Time (s)');
            title([title_figure]);
            legend(['Offset Q_y = ' num2str(struct_offsets.iyaoff)]);
            
            
            
            set(gcf,'Position',p.Results.POSITION);
            set(gcf,'PaperPosition',p.Results.PAPERPOSITION);
            
            %}
            
             % Display tauML integrated over columns:
           
        end
        
        function [] = make_subplot_log(Q,tauML,varargin)
            
            p2 = inputParser;
            addRequired(p2,'Q');
            addRequired(p2,'tauML');
            
            
            addParameter(p2,'tauML_th',[],@isnumeric);
            addParameter(p2,'iCTR',[],@isnumeric);
            addParameter(p2,'nsteps',8,@isnumeric);
            addParameter(p2,'XLABEL','Q_x (pixels)',@ischar);
            addParameter(p2,'TITLE','positive',@ischar);
            
            parse(p2,Q,tauML,varargin{:});
            
           if strcmp(p.Results.TITLE,'negative')
              Q = abs(Q); 
           end
            
             hl = line(Q , tauML);
            
            set(hl,'LineStyle','none','Marker','o','Color',[0.0 0.45 0.74]);
            pa = axis;
            
            if ~isempty(p2.Results.tauML_th)
                hl = line(Q,p2.Results.tauML_th);
                set(hl,'LineStyle','-','Color','m','LineWidth',2.0);
                axis(pa);
            end
            
            % Show CTR positions
            if ~isempty(p2.Results.iCTR)
                for ii = p2.Results.iCTR
                    hl = line(p2.Results.nsteps*ii*[1 1],pa(3:4));
                    set(hl,'LineStyle','--','Color','r');
                end
            end
               
            set(gca,'YScale','log','XScale','log');
            set(gca,'Ylim',p2.Results.YLIM);
            xlabel(p2.Results.XLABEL);
            ylabel('Mean Correlation Time (s)');
            title([p2.Results.TITLE]);
            %legend(['Offset Q_y = ' num2str(struct_offsets.iyaoff)]);
            
            
            
        end
        
        %%%% test of the SG smoothing:
        
        function [] = show_filt_Ibar(III,Ibar,dlnI,filt,maxpd,row_to_plot,frames_to_plot,figNum)
            
            %Ntb = size(III,3);
            
            
            
            figure(figNum);
            for kk = frames_to_plot%1:100:Ntb
                subplot(231);
                imagesc(log10(III(:,:,kk)));
                hold on;
                hl = line([1:size(III,2)],ones(size(III,2),1)*row_to_plot);
                set(hl,'Color','r');
                colorbar;
                title(['III at frame = ' num2str(kk)]);
                
                subplot(232);
                imagesc(filt);
                colorbar;
                title(['filter with maxpd =' num2str(maxpd)] );
                
                subplot(233);
                imagesc(log10(abs(Ibar(:,:,kk))));
                colorbar;
                title('Ibar');
                
                subplot(234);
                imagesc(dlnI(:,:,kk));
                colorbar;
                title('dlnI');
                set(gca,'Clim',[-7 7])
                
                subplot(235);
                plot(III(row_to_plot,:,kk));
                title(['III at row ' num2str(row_to_plot)]);
                
                subplot(236);
                plot(Ibar(row_to_plot,:,kk));title('Ibar');
                pause();
            end
            
            
            
        end
        
    end
end


