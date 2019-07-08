classdef DisplayTTM_plot
    % This library contains all the functions which allows to display the
    % time correlation functions
    properties(Constant)
    end
    
    
    methods(Static)
        
        % single data sets
        
       
        function [HSstruct] = display_IInormb(IIstruct,IInormb,Namefig,QvalFlag,fig_num,ImageJ,SINGLE,XCOLLabel,YROWLabel,AXISdet,INFOstr,sim_flag,CLIM)
           
            
            Nr = size(IInormb,1); 
            Nc = size(IInormb,2);
            
            
            % Read variables to plot 
            if QvalFlag == 0
                XCOLpts = IIstruct.XCOLpts; % = [1:IIstruct.Nc] - ImageJ as calculated in XCPS_read_data.TTsput_read
                YROWpts = IIstruct.YROWpts; % = [1:IIstruct.Nr] - ImageJ
                XLabelstr = XCOLLabel;
                YLabelstr = YROWLabel;
                Axis_imageFlag = 'image';
            else
                xccen = 1 + (Nc - 1)/2;
                yrcen = 1 + (Nr - 1)/2;
                Qvector_struct = XPCS_analysis.calculate_qval(xccen - ImageJ,yrcen - ImageJ,[1:Nc] - ImageJ,[1:Nr] - ImageJ,sim_flag);
                XCOLpts = Qvector_struct.del;
                YROWpts = Qvector_struct.nu;
                XLabelstr = '\Delta q\_del (1/Å)';
                YLabelstr = '\Delta q\_nu (1/Å)';
                Axis_imageFlag = 'square';
            end
            
           
            
            if ~isempty(AXISdet)
                Axislim_vect = AXISdet;
            else
                 Axislim_vect = [min(XCOLpts) max(XCOLpts) ...
                     min(YROWpts) max(YROWpts)];
            end
          
            figh = figure(fig_num);clf;
            %HSstruct.HS = pcolor(XCOLpts,YROWpts,log10(mean(IInormb,3)));
            HSstruct.HS = pcolor(XCOLpts,YROWpts,log10(IInormb(:,:,round(size(IInormb,3)/2))));
            
            axis image;
                        
            %HSstruct.HS = pcolor(XCOLpts,YROWpts,(IInormb(:,:,SINGLE+ImageJ)));
                        
            Titlestr = char(IIstruct.TITLEstuct.TITLEstr1,INFOstr);
            Namestr = [Namefig int2str(SINGLE) IIstruct.TITLEstuct.TITLEstr2];           
            Shading_mode = 'interp';
            Colorvector = CLIM;
            Colorbarflag = 1;
            flagPrettyPlot = 0;
            
            
            DisplayFunctions_XPCS.display_style(figh,'figure',Titlestr,Namestr,XLabelstr,YLabelstr,Shading_mode,Axislim_vect,flagPrettyPlot,Colorvector,Colorbarflag,Axis_imageFlag);
            
        end
        
        function [HSstruct] = display_IInormb_with_rois(Singlescan_struct,ROIS_struct,IInormb,Namefig,fig_num,ImageJ,SINGLE,XCOL,YROW,AXISdet,INFOstr,CLIM)
           
            
            % Read variables to plot
            XCOLpts = Singlescan_struct.IIstruct.XCOLpts;
            YROWpts = Singlescan_struct.IIstruct.YROWpts;
            XLabelstr = XCOL;
            YLabelstr = YROW;
            
            
            ROIS = ROIS_struct.ROIS;
            
            if ~isempty(AXISdet)
                Axislim_vect = AXISdet;
            else
                 Axislim_vect = [min(XCOLpts) max(XCOLpts) ...
                     min(YROWpts) max(YROWpts)];
            end
            
            figh = figure(fig_num);
            clf;
            
            sub1 = subplot(211);
            
            HSstruct.HS = pcolor(XCOLpts,YROWpts,log10(mean(IInormb,3)));
            [HSstruct.HROI,HSstruct.COLORORDER] = showrois(ROIS,figh,3.0);
            
                        
            Titlestr = char(Singlescan_struct.IIstruct.TITLEstuct.TITLEstr1,INFOstr);
            Namestr = [Namefig int2str(SINGLE) Singlescan_struct.IIstruct.TITLEstuct.TITLEstr2];           
            Shading_mode = 'interp';
            Colorvector = CLIM;
            Colorbarflag = 1;
            Axis_imageFlag = '';
            
            DisplayFunctions_XPCS.display_style(sub1,'subplot',Titlestr,Namestr,XLabelstr,YLabelstr,Shading_mode,Axislim_vect,0,Colorvector,Colorbarflag,Axis_imageFlag);
            
            sub2 = subplot(212);
            showrois(ROIS,figh,3.0,HSstruct.COLORORDER);
            %axis image;
            
            LEGEND = [char([ones(length(HSstruct.HROI),1)*'ROI #']) int2str([1:length(HSstruct.HROI)]')];
            legend(LEGEND);
            
            Titlestr = ['Enumeration of the ROIs and their colors'];
            Namestr = [Namefig int2str(SINGLE) ' and Enumeration of ROIs colors'];
            XLabelstr = XCOL;
            YLabelstr = YROW;
            Shading_mode = [''];
            Colorvector = [];
            Colorbarflag = 0;
            Axis_imageFlag = '';
            DisplayFunctions_XPCS.display_style(sub2,'subplot',Titlestr,Namestr,XLabelstr,YLabelstr,Shading_mode,Axislim_vect,0,Colorvector,Colorbarflag,Axis_imageFlag);
            
            
            HSstruct.DOCUclim = mDOCUclim(CLIM);
        end
        
        function [HSstruct] = show_rois_only(Read_Singlescan_struct,Singlescan_struct,fig_num,XCOLlabel,YROWlabel,AXISdet)
           
            % Read variables to plot:
            ROIS = Singlescan_struct.ROIS_struct.ROIS;
            HROI = Singlescan_struct.HS_struct.HROI;
            XCOLpts = Read_Singlescan_struct.IIstruct.XCOLpts;
            YROWpts = Read_Singlescan_struct.IIstruct.YROWpts;
            
            
            figh = figure(fig_num);
            clf;
            
            if isfield(Singlescan_struct.HS_struct,'COLORORDER')
                COLORORDER = Singlescan_struct.HS_struct.COLORORDER;
                [HSstruct.HROI,HSstruct.COLORORDER] = showrois(ROIS,gcf,3.0,COLORORDER);
            else
                [HSstruct.HROI,HSstruct.COLORORDER] = showrois(ROIS,gcf);
            end
            
            if ~isempty(AXISdet)
                Axislim_vect = AXISdet;
            else
                Axislim_vect = [min(XCOLpts) max(XCOLpts) min(YROWpts) max(YROWpts)];
            end
            
            LEGEND = [char([ones(length(HROI),1)*'ROI #']) int2str([1:length(HROI)]')];
            legend(LEGEND);
            
            Titlestr = ['Enumeration of the ROIs and their colors'];
            Namestr = ['Enumeration of ROIs colors'];
            XLabelstr = XCOLlabel;
            YLabelstr = YROWlabel;
            Shading_mode = [''];
            Colorvector = [];
            Colorbarflag = 0;
            Axis_imageFlag = 'image';

            DisplayFunctions_XPCS.display_style(figh,'figure',Titlestr,Namestr,XLabelstr,YLabelstr,Shading_mode,Axislim_vect,0,Colorvector,Colorbarflag, Axis_imageFlag);
            
        end
                      
        function [HS_sumimag_struct] =  make_summed_images(Singlescan_struct,ROIS_struct,SoXFLAG,SoYFLAG,LOGFLAG,fig_num,ImageJ,XCOLlabel,YROWlabel,AXISdet)            
            
            
            % Read variables to plot:
            IInormb = Singlescan_struct.IIstruct.IInormb;
            ROIS = ROIS_struct.ROIS;
            timeX = Singlescan_struct.IIstruct.timeX;
            lastframes = Singlescan_struct.IIstruct.lastframes;
            
            [SoY,SoX] = slicesumrois(IInormb,ROIS,ImageJ);
            
            %	Note - This does not 'normalize' the slices to the number of pixels summed
            %			However, after using function, use SoY.image{i}./SoY.norm{i}
            %				for ROI {i}
            
            figure(fig_num);
            clf;
            Numsubplots =  length(ROIS(:,1));
             
            for ii=1:length(ROIS(:,1))
                
                if SoXFLAG
                    Imagetoplot = SoX.images{ii};
                    nd = SoX.ndx{ii};
                    YROW = YROWlabel; % y label depends on the direction we are displaying
                    Name_spec = ['SoX ROIs #' int2str(ii)];
                    Title_spec = 'X';
                elseif SoYFLAG
                    Imagetoplot = SoY.images{ii};
                    nd = SoY.ndx{ii};
                    YROW = XCOLlabel;
                    Name_spec = ['SoY ROIs #' int2str(ii)];
                    Title_spec = 'Y';
                end
                
                if LOGFLAG
                    Imagetoplot = log10(Imagetoplot);
                end
                
                sub1 = subplot(ceil(sqrt(Numsubplots)),ceil(sqrt(Numsubplots)),ii);
               
                
                HS_sumimag_struct = pcolor(timeX,nd,Imagetoplot);
                makeyline(timeX(lastframes));
                
                Titlestr = char(Singlescan_struct.IIstruct.TITLEstuct.TITLEstr1(size(Singlescan_struct.IIstruct.TITLEstuct.TITLEstr1,1),:),[Singlescan_struct.IIstruct.DOCUInt, ' summed over ' Title_spec ' in ROI #' int2str(ii)]);
                XLabelstr = 'Time (s)';
                YLabelstr = YROW;
                Namestr = [Name_spec Singlescan_struct.IIstruct.TITLEstuct.TITLEstr2]
                Shading_mode = ['interp'];
                Colorvector = [];
                Colorbarflag = 1;
                Axis_imageFlag = [''];
                
                if ~isempty(AXISdet)
                    Axislim_vect = AXISdet;
                else
                    Axislim_vect = [min(timeX) max(timeX) min(nd) max(nd)];
                end
                
                DisplayFunctions_XPCS.display_style(sub1,'subplot',Titlestr,Namestr,XLabelstr,YLabelstr,Shading_mode,Axislim_vect,0,Colorvector,Colorbarflag,Axis_imageFlag);                               
            end
        end
        
        function [HS_sumpoints_struct] = display_sumpoints_1D(Read_Singlescan_struct,Single_scan,SoXFLAG,SoYFLAG,fig_num,ImageJ,XCOLlabel,YROWlabel,AXISdet)
            
          
            % Read variables to plot:
            IInormb = Read_Singlescan_struct.IIstruct.IInormb;
            ROIS = Single_scan.ROIS_struct.ROIS;
            HROI = Single_scan.HS_struct.HROI;
            COLORORDER = Single_scan.HS_struct.COLORORDER;
            [SoY,SoX] = slicesumrois(IInormb,ROIS,ImageJ);
            
            if SoXFLAG
                Imagetoplot = SoX.images;
                nd = SoX.ndx;
                YROW = YROWlabel; % y label depends on the direction we are displaying
                Name_spec = ['SoX and Sum Points '];
                Title_spec = 'XCOLS';
            elseif SoYFLAG
                Imagetoplot = SoY.images;
                nd = SoY.ndx;
                YROW = XCOLlabel;
                Name_spec = ['SoY and Sum Points '];
                Title_spec =  'YROWS';
            end
            
            figh = figure(fig_num);
            clf;
            
            HS_sumpoints_struct.HS(1) = semilogy(nd{1},sum(Imagetoplot{1}'));
            set(HS_sumpoints_struct.HS(1),'Color',COLORORDER(1,:),'LineWidth',2);
            
            
            hold on;
            for ii = 2:length(ROIS(:,1))
                HS_sumpoints_struct.HS(ii) = line(nd{ii},sum(Imagetoplot{ii}'));
                set(HS_sumpoints_struct.HS(ii),'Color',COLORORDER(ii,:),'LineWidth',2);
            end
            
            LEGEND = [char([ones(length(HROI),1)*'ROI #']) int2str([1:length(HROI)]')];
            legend(LEGEND);
            
            
            Titlestr = char(Read_Singlescan_struct.IIstruct.TITLEstuct.TITLEstr1(size(Read_Singlescan_struct.IIstruct.TITLEstuct.TITLEstr1,1),:),['summed over ' Title_spec ' in all ROI and scan pts']);
            XLabelstr = YROW;
            YLabelstr = 'Integrated int.';
            Namestr = [Name_spec Read_Singlescan_struct.IIstruct.TITLEstuct.TITLEstr2];
            Shading_mode = [''];
            Colorvector = [];
            Colorbarflag = 0;
            Axis_imageFlag = 'square';
            
            if ~isempty(AXISdet)
                Axislim_vect = AXISdet;
            else
                Axislim_vect = [ min(nd{1}) max(nd{1}) 0 max(sum(Imagetoplot{1}'))];
            end
            
            DisplayFunctions_XPCS.display_style(figh,'figure',Titlestr,Namestr,XLabelstr,YLabelstr,Shading_mode,Axislim_vect,0,Colorvector,Colorbarflag,Axis_imageFlag);
        end
                 
        function [] = sum_over_points_sets_in_SUMPOINTS(Singlescan_struct,SoXFLAG,SoYFLAG,fig_num,ImageJ,XCOLlabel,YROWlabel,AXISdet,POINTSUMS)   % sum over points sets in POINTSUMS
                        
            % Read variables to plot:
            IInormb = Singlescan_struct.IIstruct.IInormb;
            ROIS = Singlescan_struct.ROIS_struct.ROIS;
            
            [SoY,SoX] = slicesumrois(IInormb,ROIS,ImageJ);
            
            if SoXFLAG
                Imagetoplot = SoX.images;
                nd = SoX.ndx;
                XLabelstr = YROWlabel; % y label depends on the direction we are displaying
                Name_spec = ['SoX 1st ROI and Sum select points'];
                Title_spec = 'XCOLS';
            elseif SoYFLAG
                Imagetoplot = SoY.images;
                nd = SoY.ndx;
                XLabelstr = XCOLlabel;
                Name_spec = ['SoY 1st ROI and Sum selected points'];
                Title_spec = 'YROWS';
            end
            
            figure(fig_num);clf
            hold on;
            for ii=1:length(POINTSUMS(:,1))
                Ni = [POINTSUMS(ii,1): POINTSUMS(ii,2)]+ImageJ;  %use as matlab
                HL(ii) = semilogy(nd{1},sum(Imagetoplot{1}(:,Ni)'));
            end
            
            legend(addnames2matrix('sum between [', int2str(POINTSUMS),']'))
            
            Titlestr = char(Singlescan_struct.IIstruct.TITLEstuct.TITLEstr1,['summed over ' Title_spec ' in 1st ROI and between selected Spec Points']);
            YLabelstr = 'Int (arb)';
            Namestr = [Name_spec Singlescan_struct.IIstruct.TITLEstuct.TITLEstr2];
            Shading_mode = [''];
            Colorvector = [];
            Colorbarflag = 0;
            Axis_imageFlag = 'square';
            
             if ~isempty(AXISdet)
                Axislim_vect = AXISdet;
            else
                Axislim_vect = [ min(nd) max(nd) min(sum(Imagetoplot{1}')) max(sum(Imagetoplot{1}'))];
            end
            
            DisplayFunctions_XPSC.display_style(figh,'figure',Titlestr,Namestr,XLabelstr,YLabelstr,Shading_mode,Axislim_vect,Colorvector,Colorbarflag,Axis_imageFlag);
            
            
        end
        
        function [] = plot_summed_images(Singlescan_struct,ROIS_struct,IInormb,fig_num,ImageJ,CLIM,XCOLlabel,YROWlabel,AXISdet,INFOstr,POINTSUMS)
            
            %Read variables to plot
             Xsteps = Singlescan_struct.IIstruct.Xsteps;
             XCOLpts = Singlescan_struct.IIstruct.XCOLpts;
             YROWpts = Singlescan_struct.IIstruct.YROWpts;
             ROIS = ROIS_struct.ROIS;
            
            if isempty(POINTSUMS)
                POINTS=[1 length(Xsteps)]-ImageJ;
            else
                POINTS = POINTSUMS;
            end
            
            for ii = 1:length(POINTS(:,1))
                Ni = [POINTS(ii,1): POINTS(ii,2)]+ImageJ;
                NP = length(Ni);
                
                
                figh = figure(fig_num);clf;
                pcolor(XCOLpts,YROWpts,sum(IInormb(:,:,Ni),3)./(NP));
                showrois(ROIS,figh);
                axis image;
            end   
                
            Titlestr = char(Singlescan_struct.IIstruct.TITLEstuct.TITLEstr1,INFOstr);
            Namestr = ['Image summed over points/Num Points' Singlescan_struct.IIstruct.TITLEstuct.TITLEstr2];
            XLabelstr = XCOLlabel;
            YLabelstr = YROWlabel;
            Shading_mode = 'flat';
            Colorbarflag = 1;
            Colorvector = [CLIM];
            Axisimageflag = 'square';
            
            if ~isempty(AXISdet)
                Axislim_vect = AXISdet;
            else
                 Axislim_vect = [min(XCOLpts) max(XCOLpts) ...
                     min(YROWpts) max(YROWpts)];
            end
            
            DisplayFunctions_XPCS.display_style(figh,'figure',Titlestr,Namestr,XLabelstr,YLabelstr,Shading_mode,Axislim_vect,0,Colorvector,Colorbarflag,Axisimageflag);
                     
            
        end
        
        function [ROIS_struct] = plot_ROIs_as_counters(Read_Singlescan_struct,Singlescan_struct,ROIS_index,fig_num)
            
            [POSITION,PAPERPOSITION,FONTSIZE,CMAX,CLIM,XCOLlabel,YROWlabel,...
                AXISdet,DOCUclim,INFOstr] = XPCS_initialize_parameters.TTplot_parameters();
            
            
            [ImageJ,Xstepcol,SINGLE,BKG,...
                scanflag,imname,p_image,ending,POINTSUMS] = XPCS_initialize_parameters.TTsput_read_ini();
            
            % Read variables to plot:
            IInormb = Read_Singlescan_struct.IIstruct.IInormb;
            ROIS = Singlescan_struct.ROIS_struct.ROIS;
            timeX = Read_Singlescan_struct.IIstruct.timeX;
            lastframes = Read_Singlescan_struct.IIstruct.lastframes;
            COLORORDER = Singlescan_struct.HS_struct.COLORORDER;
            HROI = Singlescan_struct.HS_struct.HROI;
            
            [sumROIS,sumROISnormed] = sumrois_wNaN(IInormb,ROIS,ImageJ);

            
            figh = figure(fig_num);
            clf;     
            hold on;
            for ii = 1:size(ROIS_index,2)
                HL(ii) = plot(timeX,sumROIS(:,ROIS_index(ii))); 
                set(HL(ii),'Color',COLORORDER(ROIS_index(ii),:),'LineWidth',2);
            end
            
            
            LEGEND = [char([ones(length(HROI(ROIS_index)),1)*'ROI #']) int2str(ROIS_index')];
            legend(LEGEND);
            
            Namestr = ['ROIs as counters' Read_Singlescan_struct.IIstruct.TITLEstuct.TITLEstr2];
            XLabelstr = 'Time (s)';
            YLabelstr = 'Inorm(summed over ROI)';
            Titlestr = char(Read_Singlescan_struct.IIstruct.TITLEstuct.TITLEstr1,INFOstr,YLabelstr);
            Shading_mode = 'interp';
            Colorbarflag = 0;
            Colorvector = [];
            Axisimageflag = 'square';
            
            if ~isempty(AXISdet)
                Axislim_vect = AXISdet;
            else
                 Axislim_vect = [min(timeX) max(timeX) ...
                    0 max(max(sumROIS(:,1:length(ROIS(:,1)))))];
            end
            
            DisplayFunctions_XPCS.display_style(figh,'figure',Titlestr,Namestr,XLabelstr,YLabelstr,Shading_mode,Axislim_vect,0,Colorvector,Colorbarflag,Axisimageflag);
                    
            % update structure:
            ROIS_struct.sumROIS = sumROIS;
            ROIS_struct.sumROISnormed = sumROISnormed;
            
        end
   
        function display_style(figh,flag,Titlestr,Namestr,Xlabelstr,Ylabelstr,shadingstr,Axislim_vect,flagPrettyPlot,COLORVECTOR,Colorbarflag,axisimagestr)
            
            [POSITION,PAPERPOSITION,FONTSIZE,CMAX,CLIM,~,~,...
                ~,~,~] = XPCS_initialize_parameters.TTplot_parameters();
            
            switch flag
                case 'figure'
                    figure(figh);
                    %set(gcf,'Position',PAPERPOSITION);
                case 'subplot'
                    subplot(figh);
            end
            
            switch shadingstr
                case 'interp'
                    shading interp;
                case 'flat'
                    shading flat;
                otherwise
                    disp('Shading does not apply')
            end
            
            switch axisimagestr
                case 'image'
                    axis image;
                case 'square'
                    axis square;
                otherwise
                    disp('no specific format');
            end
            
            if ~isempty(Titlestr)
                title(Titlestr);
            end
            set(gcf,'Name',Namestr);
            xlabel(Xlabelstr);
            ylabel(Ylabelstr);
            
            set(gca,'FontSize',FONTSIZE);
            %set(gca,'Position',PAPERPOSITION);
            
            if ~isempty(COLORVECTOR)
                set(gca,'clim',COLORVECTOR);
            end
            
            set(gca,'TickDir','out');
            set(gca,'Xlim',Axislim_vect(1:2));
            set(gca,'Ylim',Axislim_vect(3:4));
            
            if Colorbarflag 
                colorbar;
            end
            colormap jet;
            
            if flagPrettyPlot
                prettyplot(gca,PAPERPOSITION);
            end
            
        end
    end
    end

