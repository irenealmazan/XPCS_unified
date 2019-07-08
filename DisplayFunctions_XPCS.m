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

