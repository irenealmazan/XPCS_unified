classdef DisplayFunctions_XPCS
    % This library contains all the functions which allows to display the
    % time correlation functions
    properties(Constant)
    end
    
    
    methods(Static)
        
        % single data sets
        
        function [] = display_singlebox_CCN2avg(ittccen,ittrcen,ittr,ittc,QvalFlag,ImageJ,fig_num)
            
             figure(fig_num);
             hold on;
            
            for icq = 1:numel(ittc)
               
                for irq = 1:numel(ittr)                 
                    
                    if QvalFlag
                        Qval_struct = XPCS_analysis.calculate_qval(ittccen,ittrcen,[ittr(1) ittr(3)+1]-ImageJ,[ittc(1) ittc(3)+1]-ImageJ);
                        
                        Xl = Qval_struct.nu(1);
                        Xh = Qval_struct.nu(2);
                        Yl = Qval_struct.del(1);
                        Yh = Qval_struct.del(2);
                    else
                        Yl = ittc(1)-ImageJ;
                        Yh = ittc(3)-ImageJ+1;
                        Xl = ittr(1)-ImageJ;
                        Xh = ittr(3)-ImageJ+1;
                    end
                    
                    
                   HL = line(([Yl Yl Yh Yh Yl]'),([Xl Xh Xh Xl Xl]'),'LineWidth',3,'Color','k');
%                     scatter(([Yl]'),([Xl ]'),3.0,'k');
%                     scatter(([Yl]'),([Xh ]'),3.0,'g');
%                     scatter(([Yh]'),([Xh ]'),3.0,'m');
%                     scatter(([Yh]'),([Xl ]'),3.0,'r');
                end
                
            end
            
        end
        
        
         function [] = display_grid_CCN2avg(ittccen,ittrcen,Ncq,Nrq,QvalFlag,iT,hwttr_allT,hwttc_allT,wrq_allT,wcq_allT,offsetcc_allT,offsetrc_allT,fig_num,ImageJ,D_ds,kvector,pixel_size,th_Bragg)
            
            % hwttr and hwttc are the half width of the box for rows and
            % columns respectively, wrw and wcq are the half number of
            % boxes per rows and conlumns and offsetcc and offsetrc are a
            % constant offset for columns and rows
            
            
            hwttr = hwttr_allT(iT);
            hwttc = hwttc_allT(iT);
            wrq = wrq_allT(iT);
            wcq = wcq_allT(iT);
            offsetcc = offsetcc_allT(iT);
            offsetrc = offsetrc_allT(iT);
            
            
             figure(fig_num);
             hold on;
            
            for icq = 1:Ncq
                offttc = (icq-wcq-1)*(2*hwttc+1);                              
               
                for irq = 1:Nrq                   
                    offttr = (irq - wrq-1)*(2*hwttr+1);
                    
                    if QvalFlag
                        Qval_struct = XPCS_analysis.calculate_qval(ittccen,ittrcen,ittccen + offttc + [-hwttc hwttc+1]+offsetcc-0.5,ittrcen + offttr + [-hwttr hwttr+1]+offsetrc - 0.5,D_ds,kvector,pixel_size,th_Bragg);
                        
                        Yl = Qval_struct.nu(1);
                        Yh = Qval_struct.nu(2);
                        Xl = Qval_struct.del(1);
                        Xh = Qval_struct.del(2);
                    else
                        Xl = ittccen + offttc +offsetcc - hwttc - ImageJ-0.5;
                        Xh = ittccen + offttc +offsetcc + hwttc+1- ImageJ-0.5;
                        Yl = ittrcen + offttr +offsetrc - hwttr- ImageJ-0.5;
                        Yh = ittrcen + offttr +offsetrc + hwttr+1- ImageJ-0.5;
                    end
                    
                    
                    HL = line(([Xl Xh Xh Xl Xl]'),([Yl Yl Yh Yh Yl]'),'LineWidth',3,'Color','k');
                    
                end
                
            end
            
         end
       
         function [] = display_mask_CCN2avg(ittccen,ittrcen,Ncq,Nrq,QvalFlag,iT,mask,fig_num,ImageJ,D_ds,kvector,pixel_size,th_Bragg)
            
            % hwttr and hwttc are the half width of the box for rows and
            % columns respectively, wrw and wcq are the half number of
            % boxes per rows and conlumns and offsetcc and offsetrc are a
            % constant offset for columns and rows
            
            
            hwttr = mask.hwttr_allT(iT);
            hwttc = mask.hwttc_allT(iT);
            wrq = mask.wrq_allT(iT);
            wcq = mask.wcq_allT (iT);
            offsetcc = mask.offsetcc_allT(iT);
            offsetrc =  mask.offsetrc_allT(iT);
            
            
            
             figure(fig_num);
             hold on;
            
            for icq = 1:Ncq
                offttc = (icq-wcq-1)*(2*hwttc+1);                              
               
                for irq = 1:Nrq                   
                    offttr = (irq - wrq-1)*(2*hwttr+1);
                    
                    if QvalFlag
                        Qval_struct = XPCS_analysis.calculate_qval(ittccen,ittrcen,ittccen + offttc + [-hwttc hwttc+1]+offsetcc-0.5,ittrcen + offttr + [-hwttr hwttr+1]+offsetrc - 0.5,D_ds,kvector,pixel_size,th_Bragg);
                        
                        Yl = Qval_struct.nu(1);
                        Yh = Qval_struct.nu(2);
                        Xl = Qval_struct.del(1);
                        Xh = Qval_struct.del(2);
                    else
                        Xl = ittccen + offttc +offsetcc - hwttc - ImageJ-0.5;
                        Xh = ittccen + offttc +offsetcc + hwttc+1- ImageJ-0.5;
                        Yl = ittrcen + offttr +offsetrc - hwttr- ImageJ-0.5;
                        Yh = ittrcen + offttr +offsetrc + hwttr+1- ImageJ-0.5;
                    end
                    
                    
                    HL = line(([Xl Xh Xh Xl Xl]'),([Yl Yl Yh Yh Yl]'),'LineWidth',3,'Color','w');
                    
                end
                
            end
            
        end
        
         
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
        
        function display_IInormbbref(IIbin_struct,boxcenterrc_struct,fig_ini,AXISdet,D_ds,kvector,pixel_size,th_Bragg)
            
           
            IInormbb = IIbin_struct.IInormbb;
            Xamountb = IIbin_struct.Xamountb;
            IInormbb_ref = IIbin_struct.IInormbb_ref;
            
            if isfield(IIbin_struct,'N_degree')
                Ndeg = IIbin_struct.N_degree;
            else
                Ndeg = 0;
            end
            
            Nr = size(IInormbb,1);
            Nc = size(IInormbb,2);
            xccen = 1 + (Nc - 1)/2;
            yrcen = 1 + (Nr - 1)/2;
            
            NumbSubplots = 1;
            counter_fig = 0;
            counter_pixel = 1;
            
            for ics = boxcenterrc_struct.offttc
                for irs = boxcenterrc_struct.offttr
                    
                     if mod(counter_pixel-1,NumbSubplots) == 0
                        fig_num = fig_ini+counter_fig;
                        fig_h = figure(fig_num);
                        clf;
                        counter_fig = counter_fig + 1;
                    end
                    
                    subh = subplot(sqrt(NumbSubplots),sqrt(NumbSubplots),counter_pixel-NumbSubplots*(counter_fig-1)); 
                    hold on;
                    plot(Xamountb',squeeze(IInormbb(irs,ics,:)),'ob');
                    plot(Xamountb',squeeze(IInormbb_ref(irs,ics,:) ),'k','LineWidth',3.0)
                    plot(Xamountb',mean(IInormbb(irs,ics,:))*ones(numel(Xamountb),1),'r','LineWidth',3.0)
            
                    % calculate corresponding qvalues            
                    Qval_struct = XPCS_analysis.calculate_qval(xccen,yrcen,ics,irs,D_ds,kvector,pixel_size,th_Bragg);

                    legend('IInormbb(ics,irs,:)',['Fit IInormbb to poly N = ' num2str(Ndeg)],'mean(IInormbb,3)');
                    
                    counter_pixel = counter_pixel  + 1;
                    
                    if ~isempty(AXISdet)
                        Axislim_vect = AXISdet;
                    else
                        Axislim_vect = [min(Xamountb) max(Xamountb) min(squeeze(IInormbb(irs,ics,:))) max(squeeze(IInormbb(irs,ics,:)))];
                    end
            
                    
                    Namestr = ['Pixel intensity vs time ' num2str(1+(counter_fig-1)*NumbSubplots) ' and ' num2str((counter_fig)*NumbSubplots) IIbin_struct.TITLEstuct.TITLEstr2];                   
                    Titlestr_1line = {'ics = ' num2str(ics) ' del = ' num2str(Qval_struct.del,'%10.3e') ' (1/A) '  ' irs = ' num2str(irs) ' nu = ' num2str(Qval_struct.nu,'%10.3e') ' (1/A)'};
                    Titlestr = {[Titlestr_1line{1:5}] [Titlestr_1line{6:end}]};
                    XLabelstr = 'Time/Xamountb (s)';
                    YLabelstr = 'Intensity';
                    Shading_mode = [''];
                    Colorvector = [];
                    Colorbarflag = 0;
                    Axisimagestr = 'square';
                    flagPrettyPlot = 0;
                    
                    DisplayFunctions_XPCS.display_style(subh,'subplot',Titlestr,Namestr,XLabelstr,YLabelstr,Shading_mode,Axislim_vect,flagPrettyPlot,Colorvector,Colorbarflag,Axisimagestr );

                    
                    
                end
            end
            
        end
      
        
        function display_IInormb_avg(IInormb_avg,indexq,flag_row_or_col,fig_ini,AXISdet)
            
           
            NumbSubplots = 1;
            counter_fig = 0;
            
            
            
            for kk = 1:numel(IInormb_avg)
                
                switch flag_row_or_col
                    
                    case 'row'
                        qtolook = [ 1:IInormb_avg(kk).Ncq_Nrq(2)];                       
                    case 'col'
                        qtolook = [ 1:IInormb_avg(kk).Ncq_Nrq(1)];                                               
                    otherwise
                        disp('please select rows or cols')
                        return;
                end
                
                for iq = qtolook
                    
                    switch flag_row_or_col
                        
                        case 'row'
                            timex = IInormb_avg(kk).scancq(indexq).scanrq(iq).timex;
                            IInormb = IInormb_avg(kk).scancq(indexq).scanrq(iq).IInormb_avg;
                            nu = IInormb_avg(kk).qvector.nu(iq);
                            del = IInormb_avg(kk).qvector.del(indexq);
                           
                        case 'col'                            
                            timex = IInormb_avg(kk).scancq(iq).scanrq(indexq).timex;
                            IInormb = IInormb_avg(kk).scancq(iq).scanrq(indexq).IInormb_avg;
                            nu = IInormb_avg(kk).qvector.nu(indexq);
                            del = IInormb_avg(kk).qvector.del(iq);
                    end
                    
                    
                    if mod(iq-1,NumbSubplots) == 0
                        fig_num = fig_ini+counter_fig;
                        fig_h = figure(fig_num);
                        clf;
                        counter_fig = counter_fig + 1;
                    end
                    
                    subh = subplot(sqrt(NumbSubplots),sqrt(NumbSubplots),iq-NumbSubplots*(counter_fig-1));            
                    plot(timex,IInormb,'LineWidth',3.0);

                    
                    if ~isempty(AXISdet)
                        Axislim_vect = AXISdet;
                    else
                        Axislim_vect = [min(timex) max(timex) min(IInormb) max(IInormb)];
                    end
            
                    
                    Namestr = ['IInormb averaged ' num2str(1+(counter_fig-1)*NumbSubplots) ' and ' num2str((counter_fig)*NumbSubplots) IInormb_avg.TITLEstruct.TITLEstr2];                   
                    Titlestr = ['q\_nu = ' num2str(nu,'%10.3e') '1/A ; q\_del = ' num2str(del,'%10.3e') '1/\AA'];
                    XLabelstr = 'Time (s)';
                    YLabelstr = 'Average II';
                    Shading_mode = ['interp'];
                    Colorvector = [];
                    Colorbarflag = 0;
                    Axisimagestr = 'square';
                    flagPrettyPlot = 0;
                    
                    DisplayFunctions_XPCS.display_style(subh,'subplot',Titlestr,Namestr,XLabelstr,YLabelstr,Shading_mode,Axislim_vect,flagPrettyPlot,Colorvector,Colorbarflag,Axisimagestr );

                    
                end 
            end         
        end
        
          
         function  display_IInormbb_inbox(IIbin_struct,IInormb_avg,iT,c_toplot,r_toplot,hwttr_allT,hwttc_allT,wrq_allT,wcq_allT,offsetcc_allT,offsetrc_allT, N_degree)
            
            hwttr = hwttr_allT(iT);
            hwttc = hwttc_allT(iT);
            wrq = wrq_allT(iT);
            wcq = wcq_allT(iT);
            offsetcc = offsetcc_allT(iT);
            offsetrc = offsetrc_allT(iT);
            
            
            Nts = size(IIbin_struct.IInormbb_ref,3);
            Ncs = size(IIbin_struct.IInormbb_ref,2);
            Nrs = size(IIbin_struct.IInormbb_ref,1);
            
            ittccen =1 + (Ncs - 1)/2 ;% index of col center
            ittrcen = 1 + (Nrs - 1)/2; % index of row center        

            
            Ncq = 2*wcq + 1;
            Nrq = 2*wrq + 1;
            
            IInormb_avg = FittingFunctions.fit_IInormb_avg_intime(IInormb_avg, N_degree);

            for icq = c_toplot%1:Ncq
                offttc = (icq-wcq-1)*(2*hwttc+1)+ offsetcc;               
                ittc = round(ittccen + offttc + [-hwttc:hwttc]) ;
               
                for irq = r_toplot%1:Nrq
                    offttr = (irq - wrq - 1)*(2*hwttr+1)+ offsetrc;
                    ittr = round(ittrcen + offttr + [-hwttr:hwttr]);
                    
                    counter = 1;
                    
                    Nsubplots = numel(ittc)*numel(ittr);
           
                    figure;
                    
                    for kkkk = ittc
                        for jjj = ittr
                            subplot(ceil(sqrt(Nsubplots)),ceil(sqrt(Nsubplots)),counter)
                            plot(IIbin_struct.Xamountb,squeeze(IIbin_struct.IInormbb(jjj,kkkk,:)./IIbin_struct.IInormbb_ref(jjj,kkkk,:)-1),'LineWidth',3.0);
                            hold on;
                            plot(IInormb_avg.scancq(icq).scanrq(irq).timex,IInormb_avg.scancq(icq).scanrq(irq).IInormb_avg./squeeze(IInormb_avg.IInormbb_avg_fit(irq,icq,:))-1,'Color','r','LineWidth',3.0);
                            %plot(IIbin_struct.Xamountb,squeeze(IIbin_struct.IInormbb_ref(jjj,kkkk,:)),'Color','r','LineWidth',3.0);
                            ylim([-2 2]);
                            counter = counter + 1;
                        end
                    end
                    
                  
                end
              
            end
           
            
         end
       
          
         function  display_IInormbb_histogram(IIbin_struct,IInormb_avg,iT,c_toplot,r_toplot,hwttr_allT,hwttc_allT,wrq_allT,wcq_allT,offsetcc_allT,offsetrc_allT, N_degree)
            
            hwttr = hwttr_allT(iT);
            hwttc = hwttc_allT(iT);
            wrq = wrq_allT(iT);
            wcq = wcq_allT(iT);
            offsetcc = offsetcc_allT(iT);
            offsetrc = offsetrc_allT(iT);
            
            
            Nts = size(IIbin_struct.IInormbb_ref,3);
            Ncs = size(IIbin_struct.IInormbb_ref,2);
            Nrs = size(IIbin_struct.IInormbb_ref,1);
            
            ittccen =1 + (Ncs - 1)/2 ;% index of col center
            ittrcen = 1 + (Nrs - 1)/2; % index of row center        

            
            Ncq = 2*wcq + 1;
            Nrq = 2*wrq + 1;
            
            IInormb_avg = FittingFunctions.fit_IInormb_avg_intime(IInormb_avg, N_degree);

            for icq = c_toplot%1:Ncq
                offttc = (icq-wcq-1)*(2*hwttc+1)+ offsetcc;               
                ittc = round(ittccen + offttc + [-hwttc:hwttc]) ;
               
                for irq = r_toplot%1:Nrq
                    offttr = (irq - wrq - 1)*(2*hwttr+1)+ offsetrc;
                    ittr = round(ittrcen + offttr + [-hwttr:hwttr]);
                    
                    counter = 1;
                    
                    Nsubplots = numel(ittc)*numel(ittr);
           
                    figure;
                    title(['Row = ' num2str(irq)]);
                    for kkkk = ittc
                        for jjj = ittr
                            delta_IInormbb(:) = squeeze(IIbin_struct.IInormbb(jjj,kkkk,:)./IIbin_struct.IInormbb_ref(jjj,kkkk,:)-1);
                            subplot(ceil(sqrt(Nsubplots)),ceil(sqrt(Nsubplots)),counter)
                            h = histogram(delta_IInormbb);
                           
                            counter = counter + 1;
                        end
                    end

                    set(gcf,'Name',['Row = ' num2str(irq)]);
                  
                end
              
            end
           
            
        end
       
         
        
        function display_CCN2avg(CCfunc,indexq,flag_row_or_col,fig_ini,AXISdet)
            
           
            NumbSubplots = 1;
            counter_fig = 0;
            
            
            
            for kk = 1:numel(CCfunc)
                
                switch flag_row_or_col
                    
                    case 'row'
                        qtolook = [ 1:CCfunc(kk).Ncq_Nrq(2)];                       
                    case 'col'
                        qtolook = [ 1:CCfunc(kk).Ncq_Nrq(1)];                                               
                    otherwise
                        disp('please select rows or cols')
                        return;
                end
                
                for iq = qtolook
                    
                    switch flag_row_or_col
                        
                        case 'row'
                            timex = CCfunc(kk).scancq(indexq).scanrq(iq).timex;
                            CCN2avg = CCfunc(kk).scancq(indexq).scanrq(iq).CCN2avg;
                            nu = CCfunc(kk).scancq(indexq).scanrq(iq).nu;
                            del = CCfunc(kk).scancq(indexq).scanrq(iq).del;
                        case 'col'                            
                            timex = CCfunc(kk).scancq(iq).scanrq(indexq).timex;
                            CCN2avg = CCfunc(kk).scancq(iq).scanrq(indexq).CCN2avg;
                            nu = CCfunc(kk).scancq(iq).scanrq(indexq).nu;
                            del = CCfunc(kk).scancq(iq).scanrq(indexq).del;
                    end
                    
                    if mod(iq-1,NumbSubplots) == 0
                        fig_num = fig_ini+counter_fig;
                        fig_h = figure(fig_num);
                        clf;
                        counter_fig = counter_fig + 1;
                    end
                    
                    subh = subplot(sqrt(NumbSubplots),sqrt(NumbSubplots),iq-NumbSubplots*(counter_fig-1));            
                    pcolor(timex,timex,CCN2avg);

                    
                    if ~isempty(AXISdet)
                        Axislim_vect = AXISdet;
                    else
                        Axislim_vect = [min(timex) max(timex) min(timex) max(timex)];
                    end
            
                    
                    Namestr = ['2times corr func between ' num2str(1+(counter_fig-1)*NumbSubplots) ' and ' num2str((counter_fig)*NumbSubplots) CCfunc.TITLEstruct.TITLEstr2];                   
                    Titlestr = ['q\_nu = ' num2str(nu,'%10.3e') '1/A ; q\_del = ' num2str(del,'%10.3e') '1/\AA'];
                    XLabelstr = 'Time (s)';
                    YLabelstr = 'Time (s)';
                    Shading_mode = ['interp'];
                    Colorvector = [];
                    Colorbarflag = 1;
                    Axisimagestr = 'square';
                    flagPrettyPlot = 0;
                    
                    DisplayFunctions_XPCS.display_style(subh,'subplot',Titlestr,Namestr,XLabelstr,YLabelstr,Shading_mode,Axislim_vect,flagPrettyPlot,Colorvector,Colorbarflag,Axisimagestr );

                    
                end 
            end         
        end
        
        function CCfunc = display_CCN2S(CCfunc,indexq,flag_row_or_col,fig_ini,plotrange)
            
          
            NumbSubplots = 1;
            counter_fig = 0;
           
            
            for kk = 1:numel(CCfunc)
                
                switch flag_row_or_col
                    
                    case 'row'
                        qtolook = [ 1:CCfunc(kk).Ncq_Nrq(2)];
                    case 'col'
                        qtolook = [ 1:CCfunc(kk).Ncq_Nrq(1)];
                    otherwise
                        disp('please select rows or cols')
                        return;
                end
                
                
                for iq = qtolook
                    
                    if mod(iq-1,NumbSubplots) == 0
                        fig_num = fig_ini+counter_fig;
                        figure(fig_num);
                        clf;
                        counter_fig = counter_fig + 1;
                    end
                    
                    switch flag_row_or_col
                        
                        case 'row'
                            time_1D = CCfunc(kk).scancq(indexq).scanrq(iq).time_1D;
                            CCNdtV = CCfunc(kk).scancq(indexq).scanrq(iq).CCNdtV;
                            nu = CCfunc(kk).scancq(indexq).scanrq(iq).nu;
                            del = CCfunc(kk).scancq(indexq).scanrq(iq).del;
                            
                            if isfield(CCfunc(kk).scancq(indexq).scanrq(iq),'CCNdtV_fit')
                                PLOTFITFLAG = 1;
                                xfit  = CCfunc(kk).scancq(indexq).scanrq(iq).CCNdtV_fit.x;
                                fitfunc = CCfunc(kk).scancq(indexq).scanrq(iq).CCNdtV_fit.fitfunc;
                                plegend = CCfunc(kk).scancq(indexq).scanrq(iq).CCNdtV_fit.plegend;
                                pout = CCfunc(kk).scancq(indexq).scanrq(iq).CCNdtV_fit.pout;
                            else

                                PLOTFITFLAG = 0;
                            end
                            
                        case 'col'                            
                            time_1D = CCfunc(kk).scancq(iq).scanrq(indexq).time_1D;
                            CCNdtV = CCfunc(kk).scancq(iq).scanrq(indexq).CCNdtV;
                            nu = CCfunc(kk).scancq(iq).scanrq(indexq).nu;
                            del = CCfunc(kk).scancq(iq).scanrq(indexq).del;
                            
                              if isfield(CCfunc(kk).scancq(iq).scanrq(indexq),'CCNdtV_fit')
                                PLOTFITFLAG = 1;
                                xfit  = CCfunc(kk).scancq(iq).scanrq(indexq).CCNdtV_fit.x;
                                fitfunc = CCfunc(kk).scancq(iq).scanrq(indexq).CCNdtV_fit.fitfunc;
                                plegend = CCfunc(kk).scancq(iq).scanrq(indexq).CCNdtV_fit.plegend;
                                pout = CCfunc(kk).scancq(iq).scanrq(indexq).CCNdtV_fit.pout;
                             else
                                PLOTFITFLAG = 0;
                             end
                    end
                    
                    
                    subh = subplot(sqrt(NumbSubplots),sqrt(NumbSubplots),iq-NumbSubplots*(counter_fig-1));
                    
                    plot(time_1D(plotrange),CCNdtV(plotrange));
                    
                    if  PLOTFITFLAG
                        
                        hold on;
                        plot(xfit,fitfunc,'r');
                        param_str = [];
                        for pp = 1:numel(plegend) param_str = [param_str ' ' plegend(pp).ptitle]; end
                        celltitle = {'nu = ' num2str(nu,'%10.3e') ' (1/A) '  ' del = ' num2str(del,'%10.3e') ' in 1/A'  param_str  ' = ' num2str(pout','%10.3e')};
                        ht = title({[celltitle{1:6}] [celltitle{7:8}] [celltitle{9:end}]});
                    else 
                    ht = title(['nu = ' num2str(nu,'%10.3e') , '  (1/A) del = ' num2str(del,'%10.3e') ' in 1/A']);
                    
                    end
                    
          
                    
                    
                    Namestr =  ['Fitted functions for irq between ' num2str(1+(counter_fig-1)*NumbSubplots) ' and ' num2str((counter_fig)*NumbSubplots) CCfunc.TITLEstruct.TITLEstr2];                   
                    Titlestr = [''];
                    XLabelstr = 'Time Delta (s)';
                    YLabelstr = 'Correlation';
                    Shading_mode = ['interp'];
                    Colorvector = [];
                    Colorbarflag = 0;
                    Axisimagestr = '';
                    flagPrettyPlot = 0;
                    Axislim_vect = [min(time_1D(plotrange)) max(time_1D(plotrange)) min(CCNdtV(plotrange)) max(CCNdtV(plotrange))];
                    DisplayFunctions_XPCS.display_style(subh,'subplot',Titlestr,Namestr,XLabelstr,YLabelstr,Shading_mode,Axislim_vect,flagPrettyPlot,Colorvector,Colorbarflag,Axisimagestr )
                    
                end
            end
        
        
      
            
        end
        
        
        function CCfunc = display_IIbinstruct_CCN2S_CCN2avg(Singlescans,indexq,flag_row_or_col,fig_ini,plotrange,AXISdet)
            
          
            NumbSubplots = 1;
            counter_fig = 0;
            
            CCfunc = Singlescans.CCN2S_struct;
            CCfunc_avg = Singlescans.CCN2avg_struct;
            IIbin_struct = Singlescans.IIbin_struct;
            
            IInormbb = IIbin_struct.IInormbb;
            Xamountb = IIbin_struct.Xamountb;
            IInormbb_ref = IIbin_struct.IInormbb_ref;
            
              
            if isfield(IIbin_struct,'N_degree')
                Ndeg = IIbin_struct.N_degree;
            else
                Ndeg = 0;
            end
            
           
            
            for kk = 1:numel(CCfunc)
                
                switch flag_row_or_col
                    
                    case 'row'
                        qtolook = [ 1:CCfunc(kk).Ncq_Nrq(2)];
                    case 'col'
                        qtolook = [ 1:CCfunc(kk).Ncq_Nrq(1)];
                    otherwise
                        disp('please select rows or cols')
                        return;
                end
                
                
                for iq = qtolook
                    
                    if mod(iq-1,NumbSubplots) == 0
                        fig_num = fig_ini+counter_fig;
                        figure(fig_num);
                        clf;
                        counter_fig = counter_fig + 1;
                    end
                    
                    switch flag_row_or_col
                        
                        case 'row'
                            irs = CCfunc_avg.boxcenterrc.offttr(iq);
                            ics = CCfunc_avg.boxcenterrc.offttc(indexq)-20;
                            
                            time_1D = CCfunc(kk).scancq(indexq).scanrq(iq).time_1D;
                            CCNdtV = CCfunc(kk).scancq(indexq).scanrq(iq).CCNdtV;
                            nu = CCfunc(kk).scancq(indexq).scanrq(iq).nu;
                            del = CCfunc(kk).scancq(indexq).scanrq(iq).del;
                            
                            if isfield(CCfunc(kk).scancq(indexq).scanrq(iq),'CCNdtV_fit')
                                PLOTFITFLAG = 1;
                                xfit  = CCfunc(kk).scancq(indexq).scanrq(iq).CCNdtV_fit.x;
                                fitfunc = CCfunc(kk).scancq(indexq).scanrq(iq).CCNdtV_fit.fitfunc;
                                plegend = CCfunc(kk).scancq(indexq).scanrq(iq).CCNdtV_fit.plegend;
                                pout = CCfunc(kk).scancq(indexq).scanrq(iq).CCNdtV_fit.pout;
                            else
                                
                                PLOTFITFLAG = 0;
                            end
                            
                            timex = CCfunc_avg(kk).scancq(indexq).scanrq(iq).timex;
                            CCN2avg = CCfunc_avg(kk).scancq(indexq).scanrq(iq).CCN2avg;
                            nu = CCfunc_avg(kk).scancq(indexq).scanrq(iq).nu;
                            del = CCfunc_avg(kk).scancq(indexq).scanrq(iq).del;
                            
                        case 'col'
                            irs = CCfunc_avg.boxcenterrc.offttr(indexq);
                            ics = CCfunc_avg.boxcenterrc.offttc(iq);
                            
                            
                            time_1D = CCfunc(kk).scancq(iq).scanrq(indexq).time_1D;
                            CCNdtV = CCfunc(kk).scancq(iq).scanrq(indexq).CCNdtV;
                            nu = CCfunc(kk).scancq(iq).scanrq(indexq).nu;
                            del = CCfunc(kk).scancq(iq).scanrq(indexq).del;
                            
                            if isfield(CCfunc(kk).scancq(iq).scanrq(indexq),'CCNdtV_fit')
                                PLOTFITFLAG = 1;
                                xfit  = CCfunc(kk).scancq(iq).scanrq(indexq).CCNdtV_fit.x;
                                fitfunc = CCfunc(kk).scancq(iq).scanrq(indexq).CCNdtV_fit.fitfunc;
                                plegend = CCfunc(kk).scancq(iq).scanrq(indexq).CCNdtV_fit.plegend;
                                pout = CCfunc(kk).scancq(iq).scanrq(indexq).CCNdtV_fit.pout;
                            else
                                PLOTFITFLAG = 0;
                            end
                            
                            
                            timex = CCfunc_avg(kk).scancq(iq).scanrq(indexq).timex;
                            CCN2avg = CCfunc_avg(kk).scancq(iq).scanrq(indexq).CCN2avg;
                            nu = CCfunc_avg(kk).scancq(iq).scanrq(indexq).nu;
                            del = CCfunc_avg(kk).scancq(iq).scanrq(indexq).del;
                    end
                    
                    subh1 = subplot(1,3,1);
                    hold on;
                    plot(Xamountb',squeeze(IInormbb(irs,ics,:)),'ob');
                    plot(Xamountb',squeeze(IInormbb_ref(irs,ics,:) ),'k','LineWidth',3.0)
                    plot(Xamountb',mean(IInormbb(irs,ics,:))*ones(numel(Xamountb),1),'r','LineWidth',3.0)
            
                    %legend('IInormbb(ics,irs,:)',['Fit IInormbb to poly N = ' num2str(Ndeg)],'mean(IInormbb,3)');
                    
                    
                    if ~isempty(AXISdet)
                        Axislim_vect = AXISdet;
                    else
                        Axislim_vect = [min(Xamountb) max(Xamountb) min(squeeze(IInormbb(irs,ics,:))) max(squeeze(IInormbb(irs,ics,:)))];
                    end
            
                    
                    Namestr = ['Pixel intensity vs time ' num2str(1+(counter_fig-1)*NumbSubplots) ' and ' num2str((counter_fig)*NumbSubplots) IIbin_struct.TITLEstuct.TITLEstr2];                   
                    %Titlestr_1line = {'ics = ' num2str(ics) ' del = ' num2str(del,'%10.3e') ' (1/A) '  ' irs = ' num2str(irs) ' nu = ' num2str(nu,'%10.3e') ' (1/A)'};
                    Titlestr_1line = {['ics = ' num2str(ics) ' irs = ' num2str(irs) ]};
                    Titlestr = Titlestr_1line;% {[Titlestr_1line(1:end)]};%{[Titlestr_1line{1:5}] [Titlestr_1line{6:end}]};
                    XLabelstr = 'Time/Xamountb (s)';
                    YLabelstr = 'Intensity';
                    Shading_mode = [''];
                    Colorvector = [];
                    Colorbarflag = 0;
                    Axisimagestr = 'square';
                    flagPrettyPlot = 0;
                    
                    DisplayFunctions_XPCS.display_style(subh1,'subplot',Titlestr,Namestr,XLabelstr,YLabelstr,Shading_mode,Axislim_vect,flagPrettyPlot,Colorvector,Colorbarflag,Axisimagestr );

                    
                    %subh = subplot(sqrt(NumbSubplots),sqrt(NumbSubplots),iq-NumbSubplots*(counter_fig-1));
                    subh2 = subplot(1,3,2);
                    pcolor(timex,timex,CCN2avg);
                    
                    
                    if ~isempty(AXISdet)
                        Axislim_vect = AXISdet;
                    else
                        Axislim_vect = [min(timex) max(timex) min(timex) max(timex)];
                    end
            
                    
                    Namestr = ['2times corr func between ' num2str(1+(counter_fig-1)*NumbSubplots) ' and ' num2str((counter_fig)*NumbSubplots) CCfunc.TITLEstruct.TITLEstr2];                   
                    Titlestr = ['q\_nu = ' num2str(nu,'%10.3e') '1/A ; q\_del = ' num2str(del,'%10.3e') '1/\AA'];
                    XLabelstr = 'Time (s)';
                    YLabelstr = 'Time (s)';
                    Shading_mode = ['interp'];
                    Colorvector = [];
                    Colorbarflag = 1;
                    Axisimagestr = 'square';
                    flagPrettyPlot = 0;
                    
                    DisplayFunctions_XPCS.display_style(subh2,'subplot',Titlestr,Namestr,XLabelstr,YLabelstr,Shading_mode,Axislim_vect,flagPrettyPlot,Colorvector,Colorbarflag,Axisimagestr );

                    
                    
                    subh3 = subplot(1,3,3);
                    plot(time_1D(plotrange),CCNdtV(plotrange),'LineWidth',3.0);
                    
                                        
                    if  PLOTFITFLAG
                        
                        hold on;
                        plot(xfit,fitfunc,'Color','r','LineWidth',3.0);
                        param_str = [];
                        for pp = 1:numel(plegend) param_str = [param_str ' ' plegend(pp).ptitle]; end
                        celltitle = {'nu = ' num2str(nu,'%10.3e') ' (1/A) '  ' del = ' num2str(del,'%10.3e') ' in 1/A'  param_str  ' = ' num2str(pout','%10.3e')};
                        ht = title({[celltitle{1:6}] [celltitle{7:8}] [celltitle{9:end}]});
                    else 
                        ht = title(['nu = ' num2str(nu,'%10.3e') , '  (1/A) del = ' num2str(del,'%10.3e') ' in 1/A']);
                    
                    end
                    
          
                    
                    
                    
                    
                    Namestr =  ['Fitted functions for irq between ' num2str(1+(counter_fig-1)*NumbSubplots) ' and ' num2str((counter_fig)*NumbSubplots) CCfunc.TITLEstruct.TITLEstr2];                   
                    Titlestr = [''];
                    XLabelstr = 'Time Delta (s)';
                    YLabelstr = 'Correlation';
                    Shading_mode = ['interp'];
                    Colorvector = [];
                    Colorbarflag = 0;
                    Axisimagestr = 'square';
                    flagPrettyPlot = 0;
                    Axislim_vect = [min(time_1D(plotrange)) max(time_1D(plotrange)) min(CCNdtV(plotrange)) max(CCNdtV(plotrange))];
                    DisplayFunctions_XPCS.display_style(subh3,'subplot',Titlestr,Namestr,XLabelstr,YLabelstr,Shading_mode,Axislim_vect,flagPrettyPlot,Colorvector,Colorbarflag,Axisimagestr )
                    
                 
                    set(gcf,'pos',[10 10 900 600]);
                end
            end
        
        
      
            
        end
        
        
        function [fitres_CCN2S] = display_fit_result(CCN2S_struct,indexq,flag_row_or_col,ppvector,figh)
            
            
           
            
            switch flag_row_or_col
                
                case 'row'
                    qtolook = [ 1:CCN2S_struct.Ncq_Nrq(2)];
                    %Namefig = ['Fit results of ' flag_row_or_col ' '  num2str(indexq) ' '  CCN2S_struct.TITLEstruct.TITLEstr2 ];
                     Namefig = ['Fit at constant del = ' num2str(CCN2S_struct.scancq(indexq).scanrq(1).del,'%10.3e') ' (1/Angstroms)  '  CCN2S_struct.TITLEstruct.TITLEstr1(size(CCN2S_struct.TITLEstruct.TITLEstr1,1),:)];
                case 'col'
                    qtolook = [ 1:CCN2S_struct.Ncq_Nrq(1)];
                   % Namefig = ['Fit results of ' flag_row_or_col ' '  num2str(indexq) ' '  CCN2S_struct.TITLEstruct.TITLEstr2 ];
                    Namefig = ['Fit at constant nu = ' num2str(CCN2S_struct.scancq(1).scanrq(indexq).nu,'%10.3e') ' (1/Angstroms)  '   CCN2S_struct.TITLEstruct.TITLEstr1(size(CCN2S_struct.TITLEstruct.TITLEstr1,1),:) ];

                otherwise
                    disp('please select rows or cols')
                    return;
            end
            
            
            for iq = qtolook 
                
                switch flag_row_or_col
                    
                    case 'row'
                        qvector(iq) = CCN2S_struct.scancq(indexq).scanrq(iq).nu;
                        CCNdtV_fit  = CCN2S_struct.scancq(indexq).scanrq(iq).CCNdtV_fit;  
                        title_tag = 'col';
                    case 'col'
                        qvector(iq) = CCN2S_struct.scancq(iq).scanrq(indexq).del;
                        CCNdtV_fit  = CCN2S_struct.scancq(iq).scanrq(indexq).CCNdtV_fit;   
                        title_tag = 'row';
                end
                
                
                for pp = 1:numel(CCNdtV_fit.pout)                  
                    fitres_CCN2S.pout(pp).val(iq) = CCNdtV_fit.pout(pp);
                    fitres_CCN2S.sigma(pp).val(iq) = CCNdtV_fit.sigp(pp);                   
                    fitres_CCN2S.plegend(pp).ptitle = CCNdtV_fit.plegend(pp).ptitle;
                end
                
                fitres_CCN2S.qvector = qvector;
                fitres_CCN2S.TITLEstruct = CCN2S_struct.TITLEstruct;
            end
            
            for pp=ppvector%1:numel(CCNdtV_fit.pout)
                
                %h = subplot(1,numel(CCNdtV_fit.pout),pp);
                fighandle = figure(figh+pp);
                errorbar(fitres_CCN2S.qvector,fitres_CCN2S.pout(pp).val,fitres_CCN2S.sigma(pp).val,'ob');
                drawnow;
                
                Namestr =  [CCN2S_struct.TITLEstruct.TITLEstr2];
                Titlestr = [char(CCNdtV_fit.plegend(pp).ptitle) ' ' title_tag ' = ' num2str(indexq)];
                
                YLabelstr = [char(CCNdtV_fit.plegend(pp).ptitle)];
                Shading_mode = [''];
                Colorvector = [];
                Colorbarflag = 0;
                Axisimagestr = 'square';
                flagPrettyPlot = 0;
                
                switch flag_row_or_col
                    
                    case 'row'
                        XLabelstr = 'qvector\_nu in 1/A';
                        switch pp
                            case 3
                                Axislim_vect = [-2e-2 2e-2 0 1.5e3];
                            otherwise
                                Axislim_vect = [-6e-3 0e-3 0 9e-3]%[-2e-2 2e-2  -abs(min(fitres_CCN2S.pout(pp).val))-1e-4 abs(max(fitres_CCN2S.pout(pp).val))+1e-4];
                        end
                        
                    case 'col'
                        XLabelstr = 'qvector\_del in 1/A';
                        switch pp
                            case 3
                                Axislim_vect = [-5e-4 5e-4 0 2.5e3];                                
                            otherwise
                                Axislim_vect = [-5e-4 5e-4  -abs(min(fitres_CCN2S.pout(pp).val))-1e-4 abs(max(fitres_CCN2S.pout(pp).val))+1e-4];
                        end
                end
                
                DisplayFunctions_XPCS.display_style(fighandle,'figure',Titlestr,Namestr,XLabelstr,YLabelstr,Shading_mode,Axislim_vect,flagPrettyPlot,Colorvector,Colorbarflag,Axisimagestr )
                
            end
            
            %set(gcf,'Name',Namefig);
            
%             
        end
        
        
        function display_fit_result_allcol_row(pout_struct,flag_row_or_col,ppvector,clrs,figh)
            
            % read fit results:
            for pp = ppvector
                
                for iq = [1:numel( pout_struct) ]
                    fighandle = figure(figh+pp+iq);
                    clf;
                    errorbar(pout_struct(iq).qvector,pout_struct(iq).pout(:,pp),pout_struct(iq).sigma(:,pp),clrs(2*(iq-1)+1:2*(iq-1)+2));
                    drawnow;
                    
                    
                    Namestr =  [pout_struct(iq).pout_struct.TITLEstruct.TITLEstr2];
                    
                    YLabelstr = [char(pout_struct(iq).plegend(pp).ptitle)];
                    Shading_mode = [''];
                    Colorvector = [];
                    Colorbarflag = 0;
                    Axisimagestr = 'square';
                    flagPrettyPlot = 0;
                    
                    switch flag_row_or_col
                        
                        case 'row'
                            XLabelstr = 'qvector\_nu in 1/A';
                            switch pp
                                case 3
                                    Axislim_vect = [-1e-2 1e-2 0 2.5e3];
                                otherwise
                                    % Axislim_vect = [-2e-2 2e-2  -abs(min( pout_struct.pout(:,pp)))-1e-4 abs(max( pout_struct.pout(:,pp)))+1e-4];
                                    Axislim_vect = [-5e-4 5e-4 -1 1];
                            end
                            titletag = 'col';
                            
                        case 'col'
                            XLabelstr = 'qvector\_del in 1/A';
                            switch pp
                                case 3
                                    Axislim_vect = [-5e-4 5e-4 0 2.5e3];
                                otherwise
                                    % Axislim_vect = [-5e-4 5e-4  -abs(min( pout_struct.pout(:,pp)))-1e-4 abs(max( pout_struct.pout(:,pp)))+1e-4];
                                    Axislim_vect = [-5e-4 5e-4 -1 1];
                            end
                            titletag = 'row';
                    end
                    
                    Titlestr = [char(pout_struct(iq).plegend(pp).ptitle) '  ' titletag ' = ' num2str(iq)];

                    DisplayFunctions_XPCS.display_style(fighandle,'figure',Titlestr,Namestr,XLabelstr,YLabelstr,Shading_mode,Axislim_vect,flagPrettyPlot,Colorvector,Colorbarflag,Axisimagestr )
                    
                    
                end
                
                
                            
            end
            
        end
        
        function display_fit_result_allcol_row_inonefig(pout_struct,flag_row_or_col,ppvector,clrs,figh)
            
            % read fit results:
            for pp = ppvector
                fighandle = figure(figh+pp);
                    clf;
                    hold on;
                for iq = [1:numel( pout_struct) ]
                    
                    errorbar(pout_struct(iq).qvector,pout_struct(iq).pout(:,pp),pout_struct(iq).sigma(:,pp),['-' clrs(2*(iq-1)+1:2*(iq-1)+2)]);
                    drawnow;
                    
                    
                    Namestr =  [pout_struct(iq).TITLEstruct.TITLEstr2];
                    
                    YLabelstr = [char(pout_struct(iq).plegend(pp).ptitle)];
                    Shading_mode = [''];
                    Colorvector = [];
                    Colorbarflag = 0;
                    Axisimagestr = 'square';
                    flagPrettyPlot = 0;
                    
                    switch flag_row_or_col
                        
                        case 'row'
                            XLabelstr = 'qvector\_nu in 1/A';
                            switch pp
                                case 3
                                    Axislim_vect = [-1e-2 1e-2 0 2.5e3];
                                otherwise
                                    % Axislim_vect = [-2e-2 2e-2  -abs(min( pout_struct.pout(:,pp)))-1e-4 abs(max( pout_struct.pout(:,pp)))+1e-4];
                                    Axislim_vect = [-5e-4 5e-4 -1 1];
                            end
                            titletag = 'col';
                            
                        case 'col'
                            XLabelstr = 'qvector\_del in 1/A';
                            switch pp
                                case 3
                                    Axislim_vect = [-5e-4 5e-4 0 2.5e3];
                                otherwise
                                    % Axislim_vect = [-5e-4 5e-4  -abs(min( pout_struct.pout(:,pp)))-1e-4 abs(max( pout_struct.pout(:,pp)))+1e-4];
                                    Axislim_vect = [-5e-4 5e-4 -1 1];
                            end
                            titletag = 'row';
                    end
                    
                    Titlestr = [char(pout_struct(iq).pout_struct.plegend(pp).ptitle) '  ' titletag ' = ' num2str(iq)];

                    DisplayFunctions_XPCS.display_style(fighandle,'figure',Titlestr,Namestr,XLabelstr,YLabelstr,Shading_mode,Axislim_vect,flagPrettyPlot,Colorvector,Colorbarflag,Axisimagestr )
                    
                    
                end
                
                
                            
            end
            
        end
        
        function pout_surf_struct = display_fit_result_allcol_row_pcolor(Allscans,flag_row_or_col,ppvector,figh)
            
            pout_struct = Allscans.fitres_CCN2S;
            qvector_struct = Allscans.qvector_struct;
            % read fit results:
            for pp = ppvector
                fighandle = figure(figh+pp);
                clf;
                hold on;
                
              
                for iq = [1:numel( pout_struct) ]
                    
                   
                    
                    switch flag_row_or_col
                        case 'col'
                            pout_surf_struct(pp).pout(iq,:) =  pout_struct(iq).pout(pp).val;
                            pout_surf_struct(pp).sigma(iq,:) =  pout_struct(iq).sigma(pp).val;
                            
                        case 'row'
                            pout_surf_struct(pp).pout(:,iq) =   pout_struct(iq).pout(pp).val;
                            pout_surf_struct(pp).sigma(:,iq) =   pout_struct(iq).sigma(pp).val;
                    end
                    
                end
                pout_surf_struct(pp).del = [qvector_struct.del];
                pout_surf_struct(pp).nu = [qvector_struct.nu];
                pout_surf_struct(pp).ptitle = pout_struct(iq).plegend(pp).ptitle;
                pout_surf_struct(pp).TITLEstruct = pout_struct(iq).TITLEstruct;
                
                subh1 = subplot(1,2,1);
                %pcolor(pout_surf_struct(pp).del,pout_surf_struct(pp).nu,pout_surf_struct(pp).pout);
                imagesc(pout_surf_struct(pp).del,pout_surf_struct(pp).nu,pout_surf_struct(pp).pout);
                
                %%{
                Titlestr = [char(pout_surf_struct(pp).ptitle)];
                Namestr = '';
                XLabelstr = 'qvector\_del in 1/A';    
                YLabelstr =  'qvector\_nu in 1/A';   
                Shading_mode = [''];
                
                if pp == 3
                    Colorvector = [0 3e3];
                else
                    Colorvector = [-3e-2 3e-2];
                end
                
                Colorbarflag = 1;
                Axisimagestr = 'square';
                flagPrettyPlot = 0;
                Axislim_vect = [-6e-4 6e-4 -1e-2 1e-2];
                
                DisplayFunctions_XPCS.display_style(subh1,'subplot',Titlestr,Namestr,XLabelstr,YLabelstr,Shading_mode,Axislim_vect,flagPrettyPlot,Colorvector,Colorbarflag,Axisimagestr )
%}
                
                %%{
                subh2 = subplot(1,2,2);
                %index_nan = isnan(real(pout_surf_struct(pp).sigma));
                %index_no_nan = find(index_nan == 0);
                 imagesc(pout_surf_struct(pp).del,pout_surf_struct(pp).nu,real(pout_surf_struct(pp).sigma));
                
                Titlestr = [char(pout_surf_struct(pp).ptitle) 'sigma'];
                Namestr = '';
                XLabelstr = 'qvector\_del in 1/A';    
                YLabelstr =  'qvector\_nu in 1/A';   
                Shading_mode = [''];
                if pp == 3
                    Colorvector = [0 3e3];
                else
                     Colorvector = [-3e-2 3e-2];
                end
                Colorbarflag = 1;
                Axisimagestr = 'square';
                flagPrettyPlot = 0;
                Axislim_vect = [-6e-4 6e-4 -1e-2 1e-2];
                
                DisplayFunctions_XPCS.display_style(subh2,'subplot',Titlestr,Namestr,XLabelstr,YLabelstr,Shading_mode,Axislim_vect,flagPrettyPlot,Colorvector,Colorbarflag,Axisimagestr )
%}
                
                
                Namestr =  [pout_struct(1).TITLEstruct.TITLEstr2];
                set(gcf,'Name',Namestr);
            end
           
        end
       
        
       function [] = display_fit_average_results(pout_avg,ppvector,figh)
            
          
            
            % read fit results:
            for pp = ppvector
                fighandle = figure(figh+pp);
                clf;
                hold on;
              
              
                subh1 = subplot(1,2,1);
                %axes('box','on');
                plot(pout_avg(pp).del,pout_avg(pp).weight_mean_row);
                title(['parameter = ' num2str(pp) ' row mean']);
                xlabel('del (1/Angstroms)');
                
                subh2 = subplot(1,2,2);
                %axes('box','on');
                plot(pout_avg(pp).nu,pout_avg(pp).weight_mean_col);
                title(['parameter = ' num2str(pp) ' col. mean']);
                xlabel('nu (1/Angstroms)');
                
                %subh1 = subplot(1,2,1);
                %pcolor(pout_surf_struct(pp).del,pout_surf_struct(pp).nu,pout_surf_struct(pp).pout);
                %{
                Titlestr = [char(pout_struct(1).pout_struct.plegend(pp).ptitle)];
                Namestr = '';
                XLabelstr = 'qvector\_del in 1/A';    
                YLabelstr =  'qvector\_nu in 1/A';   
                Shading_mode = [''];
                
                if pp == 3
                    Colorvector = [0 3e3];
                else
                    Colorvector = [-3e-2 3e-2];
                end
                
                Colorbarflag = 1;
                Axisimagestr = 'square';
                flagPrettyPlot = 0;
                Axislim_vect = [-6e-4 6e-4 -1e-2 1e-2];
                
                DisplayFunctions_XPCS.display_style(subh1,'subplot',Titlestr,Namestr,XLabelstr,YLabelstr,Shading_mode,Axislim_vect,flagPrettyPlot,Colorvector,Colorbarflag,Axisimagestr )
%}
                
                %{
                subh2 = subplot(1,2,2);
                pcolor(pout_surf_struct(pp).del,pout_surf_struct(pp).nu,real(pout_surf_struct(pp).sigma));
                
                Titlestr = [char(pout_struct(1).pout_struct.plegend(pp).ptitle) 'sigma'];
                Namestr = '';
                XLabelstr = 'qvector\_del in 1/A';    
                YLabelstr =  'qvector\_nu in 1/A';   
                Shading_mode = [''];
                if pp == 3
                    Colorvector = [0 3e3];
                else
                     Colorvector = [-3e-2 3e-2];
                end
                Colorbarflag = 1;
                Axisimagestr = 'square';
                flagPrettyPlot = 0;
                Axislim_vect = [-6e-4 6e-4 -1e-2 1e-2];
                
                DisplayFunctions_XPCS.display_style(subh2,'subplot',Titlestr,Namestr,XLabelstr,YLabelstr,Shading_mode,Axislim_vect,flagPrettyPlot,Colorvector,Colorbarflag,Axisimagestr )
%}
                
                
                Namestr =  [pout_avg(pp).TITLEstruct.TITLEstr2];
                set(gcf,'Name',Namestr);
            end
           
        end
       
     
       function pout_surf_struct = display_fit_result_allcol_row_pcolor_weight(Allscan_fit,flag_row_or_col,ppvector,figh)
            
            pout_struct = Allscan_fit.pout_all;
            qvector_struct = Allscan_fit.qvector_struct;
            
            % read fit results:
            for pp = ppvector
                fighandle = figure(figh+pp);
                clf;
                hold on;
                
                for iq = [1:numel( pout_struct) ]
                    
                    switch flag_row_or_col
                        case 'col'
                            pout_surf_struct(pp).pout(iq,:) =  pout_struct(iq).pout_struct.pout(:,pp);
                            pout_surf_struct(pp).sigma(iq,:) =  pout_struct(iq).pout_struct.sigma(:,pp);
                            
                        case 'row'
                            pout_surf_struct(pp).pout(:,iq) =  pout_struct(iq).pout_struct.pout(:,pp);
                            pout_surf_struct(pp).sigma(:,iq) =  pout_struct(iq).pout_struct.sigma(:,pp);
                    end
                    
                end
                pout_surf_struct(pp).nu =qvector_struct.nu;
                pout_surf_struct(pp).del =qvector_struct.del;
                 
                pcolor(pout_surf_struct(pp).del,pout_surf_struct(pp).nu,pout_surf_struct(pp).pout./real(pout_surf_struct(pp).sigma));
                
                Titlestr = [char(pout_struct(1).pout_struct.plegend(pp).ptitle)];
                Namestr = '';
                XLabelstr = 'qvector\_del in 1/A';    
                YLabelstr =  'qvector\_nu in 1/A';   
                Shading_mode = [''];
                
                if pp == 3
                    Colorvector = [0 5e1];
                else
                    Colorvector = [-3e-2 3e-2];
                end
                
                Colorbarflag = 1;
                Axisimagestr = 'square';
                flagPrettyPlot = 0;
                Axislim_vect = [-6e-4 6e-4 -1e-2 1e-2];                               
                Namestr =  [pout_struct(1).pout_struct.TITLEstruct.TITLEstr2];
                
                DisplayFunctions_XPCS.display_style(fighandle,'figure',Titlestr,Namestr,XLabelstr,YLabelstr,Shading_mode,Axislim_vect,flagPrettyPlot,Colorvector,Colorbarflag,Axisimagestr )
            end
           
        end
        
       
       function [pout_struct,fig_handle] = display_fit_result_log(Allscans,pp_index,indexq,flag_row_or_col,figh)
            
            fig_handle = figure(figh);
            clf;
            
            NumofSubplots = numel(Allscans);
            
            for iT = 1:numel(Allscans)
                
                  switch flag_row_or_col
                    
                    case 'row'
                        qtolook = [ 1:Allscans(iT).CCN2S_struct.Ncq_Nrq(2)];
                    case 'col'
                        qtolook = [ 1:Allscans(iT).CCN2S_struct.Ncq_Nrq(1)];
                    otherwise
                        disp('please select rows or cols')
                        return;
                end
                
                
                counter_neg = 1;
                counter_pos = 1;

                %subh = subplot(1,NumofSubplots,iT);
                hold on;
                
                for iq = qtolook
                    
                    switch flag_row_or_col
                        
                        case 'row'
                            qvector(iq) = Allscans(iT).CCN2S_struct.scancq(indexq).scanrq(iq).nu;
                            pout  = Allscans(iT).CCN2S_struct.scancq(indexq).scanrq(iq).CCNdtV_fit.pout(pp_index);
                            sigp  = Allscans(iT).CCN2S_struct.scancq(indexq).scanrq(iq).CCNdtV_fit.sigp(pp_index);
                            ptitle = Allscans(iT).CCN2S_struct.scancq(indexq).scanrq(iq).CCNdtV_fit.plegend(pp_index).ptitle;
                            Axislim_vect = [3e-5 3e-2 1 1e4];
                            XLabelstr = 'qvector\_nu in 1/A';
                        case 'col'
                            qvector(iq) = Allscans(iT).CCN2S_struct.scancq(iq).scanrq(indexq).del;
                            pout  = Allscans(iT).CCN2S_struct.scancq(iq).scanrq(indexq).CCNdtV_fit.pout(pp_index);
                            sigp  = Allscans(iT).CCN2S_struct.scancq(iq).scanrq(indexq).CCNdtV_fit.sigp(pp_index);
                            ptitle = Allscans(iT).CCN2S_struct.scancq(iq).scanrq(indexq).CCNdtV_fit.plegend(pp_index).ptitle;
                            Axislim_vect = [1e-5 1e-3 1 1e4];
                            XLabelstr = 'qvector\_del in 1/A';
                    end
                    
                    
                    
                    if qvector(iq)== 0
                        qvector(iq) = -1e-6;
                    end
                    
                    
                    if qvector(iq)< 0
                        pout_negative(counter_neg) = pout;
                        sigma_negative(counter_neg) = sigp;
                        qvector_toplot_negative(counter_neg) = - qvector(iq);
                        
                        
                        counter_neg = counter_neg + 1;
                    else
                        pout_positive(counter_pos) = pout;
                        sigma_positive(counter_pos) = sigp;
                        qvector_toplot_positive(counter_pos) = qvector(iq);
                        
                        counter_pos = counter_pos + 1;
                    end
                end
                
                
                if counter_neg == 1
                    pout_struct(iT).tau_negative = ones(1,1);%[pout_negative pout_positive];
                    pout_struct(iT).sigma_negative = ones(1,1);%[sigma_negative sigma_positive];
                    pout_struct(iT).qvector_negative = ones(1,1);%[qvector_toplot_negative qvector_topl
                else
                    pout_struct(iT).tau_negative = [pout_negative];%[pout_negative pout_positive];
                    pout_struct(iT).sigma_negative = [sigma_negative];%[sigma_negative sigma_positive];
                    pout_struct(iT).qvector_negative = [qvector_toplot_negative];%[qvector_toplot_negative qvector_toplot_positive];
                    pout_struct(iT).fitrange_negative = find(~isnan(pout_struct(iT).sigma_negative));
                    qvector_compatible_negative = pout_struct(iT).qvector_negative(find(~isnan(pout_struct(iT).sigma_negative)));
                    pout_struct(iT).fitrange_negative = pout_struct(iT).fitrange_negative(find(qvector_compatible_negative~=0));
                end
                if counter_pos == 1
                    pout_struct(iT).tau_positive = ones(1,1);%[pout_negative pout_positive];
                    pout_struct(iT).sigma_positive = ones(1,1);
                    pout_struct(iT).qvector_positive = ones(1,1);
                else
                    pout_struct(iT).tau_positive = [pout_positive];%[pout_negative pout_positive];
                    pout_struct(iT).sigma_positive = [sigma_positive];%[sigma_negative sigma_positive];
                    pout_struct(iT).qvector_positive = [qvector_toplot_positive];%[qvector_toplot_negative qvector_toplot_positive];
                    pout_struct(iT).fitrange_positive = find(~isnan(pout_struct(iT).sigma_positive));
                    qvector_compatible_positive = pout_struct(iT).qvector_positive(find(~isnan(pout_struct(iT).sigma_positive)));
                    pout_struct(iT).fitrange_positive = pout_struct(iT).fitrange_positive(find(qvector_compatible_positive~=0));
                end
                pout_struct(iT).plegend(pp_index).ptitle = ptitle;
                pout_struct(iT).TITLEstruct = Allscans(iT).CCN2S_struct.TITLEstruct;
                
                subh1 = subplot(212);
                errorbar(pout_struct(iT).qvector_negative,pout_struct(iT).tau_negative,pout_struct(iT).sigma_negative,'*k');
                %errorbar(qvector_toplot_negative,pout_negative,sigma_negative,'ob');
                %errorbar(qvector_toplot_positive,pout_positive,sigma_positive,'or');
                set(gca,'Yscale','log');
                set(gca,'Xscale','log');
                
                
                
                 Namestr =  [ptitle];
                
                switch flag_row_or_col
                    
                    case 'row'
                        Titlestr = {[Allscans(iT).CCN2S_struct.TITLEstruct.TITLEstr1(size(Allscans(iT).CCN2S_struct.TITLEstruct.TITLEstr1,1),:)  ['col = ' num2str(indexq) 'negative']]};
                    case 'col'
                        Titlestr = {[Allscans(iT).CCN2S_struct.TITLEstruct.TITLEstr1(size(Allscans(iT).CCN2S_struct.TITLEstruct.TITLEstr1,1),:)  ['row = ' num2str(indexq) 'negative']]};
                end
            
                
                YLabelstr = ptitle;
                Shading_mode = ['interp'];
                Colorvector = [];
                Colorbarflag = 0;
                Axisimagestr = 'square';
                flagPrettyPlot = 0;
               
                
                DisplayFunctions_XPCS.display_style(subh1,'subplot',Titlestr,Namestr,XLabelstr,YLabelstr,Shading_mode,Axislim_vect,flagPrettyPlot,Colorvector,Colorbarflag,Axisimagestr )
                
                
                subh2 = subplot(211);
                errorbar(pout_struct(iT).qvector_positive,pout_struct(iT).tau_positive,pout_struct(iT).sigma_positive,'*b');
                %errorbar(qvector_toplot_negative,pout_negative,sigma_negative,'ob');
                %errorbar(qvector_toplot_positive,pout_positive,sigma_positive,'or');
                set(gca,'Yscale','log');
                set(gca,'Xscale','log');
                
                %legend('negative q','positive q');
                
                Namestr =  [ptitle];
                
                switch flag_row_or_col
                    
                    case 'row'
                        Titlestr = {[Allscans(iT).CCN2S_struct.TITLEstruct.TITLEstr1(size(Allscans(iT).CCN2S_struct.TITLEstruct.TITLEstr1,1),:)  ['col = ' num2str(indexq) 'positive']]};
                    case 'col'
                        Titlestr = {[Allscans(iT).CCN2S_struct.TITLEstruct.TITLEstr1(size(Allscans(iT).CCN2S_struct.TITLEstruct.TITLEstr1,1),:)  ['row = ' num2str(indexq) 'positive']]};
                end
            
                
                YLabelstr = ptitle;
                Shading_mode = ['interp'];
                Colorvector = [];
                Colorbarflag = 0;
                Axisimagestr = 'square';
                flagPrettyPlot = 0;
               
                
                DisplayFunctions_XPCS.display_style(subh2,'subplot',Titlestr,Namestr,XLabelstr,YLabelstr,Shading_mode,Axislim_vect,flagPrettyPlot,Colorvector,Colorbarflag,Axisimagestr )
                
              
                
            end
        end
            
       
       function [pout_struct] = append_q_law(figh,pout_struct,SplitQ_CTR,t_osc,color)
            
            slope = SplitQ_CTR*t_osc/(2*pi); % see formula 8 in the overleaf document about step_velocity_XPCS
            
            
            for iT = 1:numel(pout_struct)
                 qvector_negative = pout_struct(iT).qvector_negative;
                 tau_theory_negative = slope./qvector_negative;
                 pout_struct(iT).tau_theory_negative = tau_theory_negative;
                 
                 qvector_positive = pout_struct(iT).qvector_positive;
                 tau_theory_positive = slope./qvector_positive;
                 pout_struct(iT).tau_theory_negative = tau_theory_positive;
                 pout_struct(iT).slope =slope;
            end
            
            figure(figh);
            
            subplot(212);
            hold on;
            plot(qvector_negative,tau_theory_negative,'Color',color,'LineWidth',3.0);
            
            subplot(211);
            hold on;
            plot(qvector_positive,tau_theory_positive,'Color',color,'LineWidth',3.0);
            
        end
        
        % all data sets:
        
       function [HSstruct] = display_IInormb_all(Allscans,Namefig,QvalFlag,LOGFLAG,fig_ini)
            
            % Read predefined parameters:
            [ImageJ,~,SINGLE,~,~,~,...
                ~,~,~] = XPCS_initialize_parameters.TTsput_read_ini();
            
            [~,~,~,~,CLIM,XCOL,...
                YROW,AXISdet,~,INFOstr] = XPCS_initialize_parameters.TTplot_parameters();
            
            NumofSubplots = 3;
            NumofFigs = round(numel(Allscans)/NumofSubplots);
            counter_fig = 0;
            
            for iT = 1:numel(Allscans)
                
                if ~LOGFLAG
                    IInormb = Allscans(iT).IIstruct.IInormb;
                else
                    IInormb = log10(Allscans(iT).IIstruct.IInormb);
                    Namefig = [Namefig 'in log'];
                end
                
                Nr = size(IInormb,2);
                Nc = size(IInormb,1);
                
                
                % Read variables to plot
                if QvalFlag == 0
                    XCOLpts = Allscans(iT).IIstruct.XCOLpts;
                    YROWpts = Allscans(iT).IIstruct.YROWpts;
                    XLabelstr = XCOL;
                    YLabelstr = YROW;
                    Axis_imageFlag = 'image';
                else
                    xccen = 1 + (Nc - 1)/2;
                    yrcen = 1 + (Nr - 1)/2;
                    Qvector_struct = XPCS_analysis.calculate_qval(xccen,yrcen,[1:Nc],[1:Nr]);
                    XCOLpts = Qvector_struct.del;
                    YROWpts = Qvector_struct.nu;
                    XLabelstr = '\Delta q\_del';
                    YLabelstr = '\Delta q\_nu';
                    Axis_imageFlag = 'square';
                end
                
                
                
                if ~isempty(AXISdet)
                    Axislim_vect = AXISdet;
                else
                    Axislim_vect = [min(XCOLpts) max(XCOLpts) 70 160];
                end
                
                if mod(iT-1,NumofSubplots) == 0
                    fig_num = fig_ini+counter_fig;
                    figure(fig_num);
                    clf;
                    set(gcf,'Name',Namefig);
                    counter_fig = counter_fig + 1;
                end
                
                subh = subplot(NumofSubplots,1,iT-NumofSubplots*(counter_fig-1));
                HSstruct.HS = pcolor(XCOLpts,YROWpts,(IInormb(:,:,SINGLE+ImageJ)));
                
                Titlestr = char(Allscans(iT).IIstruct.TITLEstuct.TITLEstr1(3,23:26),INFOstr);
                Namestr = [''];
                Shading_mode = 'interp';
                Colorvector = [];
                Colorbarflag = 1;
                flagPrettyPlot = 0;
                
                
                DisplayFunctions_XPCS.display_style(subh,'subplot',Titlestr,Namestr,XLabelstr,YLabelstr,Shading_mode,Axislim_vect,flagPrettyPlot,Colorvector,Colorbarflag,Axis_imageFlag);
                
                
                
            end
            
            
            
           
        end
        
       function vel_struct = display_vel_vs_temp(TCV,pout_struct,Namefig,island_size,fignum)
            
            
            for iT = 1:numel(pout_struct)
                vel_struct.temp(iT) = TCV(iT);
                
                vel_struct.vel_AngsperSec(iT) = 1/pout_struct(iT).pout_fit.pout.pout;
                vel_struct.vel_sigma_AngsperSec(iT) = abs(pout_struct(iT).pout_fit.pout.sigp/pout_struct(iT).pout_fit.pout.pout^2);
                vel_struct.growth_invSec(iT) = vel_struct.vel_AngsperSec(iT)/island_size;
                vel_struct.growth_sigma_invSec(iT) = vel_struct.vel_sigma_AngsperSec(iT)/island_size;
            end
            
            fig_h = figure(fignum);            
            subh1 = subplot(1,2,1);
            errorbar(vel_struct.temp,vel_struct.vel_AngsperSec,vel_struct.vel_sigma_AngsperSec,'ob');
            subh2 = subplot(1,2,2);
            errorbar(vel_struct.temp,vel_struct.growth_invSec,vel_struct.growth_sigma_invSec,'ob');
            
            
            
            
            Namestr =  [''];
            Titlestr = ['v vs temperature'];
            XLabelstr = 'Temperature (C)';
            YLabelstr = 'v (Angstroms/s)';
            Shading_mode = [''];
            Colorvector = [];
            Colorbarflag = 0;
            Axisimagestr = 'square';
            flagPrettyPlot = 0;
            Axislim_vect = [vel_struct.temp(1)-100 vel_struct.temp(end)+100 0 10];
            
            
            DisplayFunctions_XPCS.display_style(subh1 ,'subplot',Titlestr,Namestr,XLabelstr,YLabelstr,Shading_mode,Axislim_vect,flagPrettyPlot,Colorvector,Colorbarflag,Axisimagestr )
            
            Namestr =  [Namefig];
            Titlestr = ['Growth rate vs temperature'];
            XLabelstr = 'Temperature (C)';
            YLabelstr = 'Growth rate (1/sec)';
            Shading_mode = [''];
            Colorvector = [];
            Colorbarflag = 0;
            Axisimagestr = 'square';
            flagPrettyPlot = 0;
            Axislim_vect = [vel_struct.temp(1)-100 vel_struct.temp(end)+100  0 4e-3];
            
            
            DisplayFunctions_XPCS.display_style(subh2 ,'subplot',Titlestr,Namestr,XLabelstr,YLabelstr,Shading_mode,Axislim_vect,flagPrettyPlot,Colorvector,Colorbarflag,Axisimagestr )
            
            
            
            

        end
        
         function display_growth_vs_power(power_vector,vel_struct,index_vect, temperature_str,fignum)
            
            figure(fignum);
            subh1 = subplot(1,2,1);
            errorbar(power_vector,vel_struct.vel_AngsperSec(index_vect),vel_struct.vel_sigma_AngsperSec(index_vect),'ob');
            
            subh2 = subplot(1,2,2);
            errorbar(power_vector,vel_struct.growth_invSec(index_vect),vel_struct.growth_sigma_invSec(index_vect),'ob');
            
            Namestr =  ['Growth rate and velocity vs power at ' temperature_str];
            Titlestr = ['Velocity vs power at ' temperature_str];
            XLabelstr = 'Power W';
            YLabelstr = 'Velocity (A/sec)';
            Shading_mode = [''];
            Colorvector = [];
            Colorbarflag = 0;
            Axisimagestr = 'square';
            flagPrettyPlot = 0;
            Axislim_vect = [power_vector(1)-2 power_vector(end)+2 0 20];
            
            
            DisplayFunctions_XPCS.display_style( subh1 ,'subplot',Titlestr,Namestr,XLabelstr,YLabelstr,Shading_mode,Axislim_vect,flagPrettyPlot,Colorvector,Colorbarflag,Axisimagestr )
            
            
            
            
            Namestr =  ['Growth rate and velocity vs power at ' temperature_str];
            Titlestr = ['Growth rate vs power at ' temperature_str];
            XLabelstr = 'Power W';
            YLabelstr = 'Growth rate (1/sec)';
            Shading_mode = [''];
            Colorvector = [];
            Colorbarflag = 0;
            Axisimagestr = 'square';
            flagPrettyPlot = 0;
            Axislim_vect = [power_vector(1)-2 power_vector(end)+2 0 1e-2];
            
            
            DisplayFunctions_XPCS.display_style( subh2 ,'subplot',Titlestr,Namestr,XLabelstr,YLabelstr,Shading_mode,Axislim_vect,flagPrettyPlot,Colorvector,Colorbarflag,Axisimagestr )
            

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

