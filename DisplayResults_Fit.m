classdef DisplayResults_Fit
    % This library contains all the functions which allows to display the
    % time correlation functions
    properties(Constant)
    end
    
    
    methods(Static)
        
        
        function [fitres_CCN2S] = display_fit_result(CCN2S_struct,indexq,flag_row_or_col,ppvector,figh)
            
            switch flag_row_or_col
                
                case 'row'
                    qtolook = [ 1:CCN2S_struct.Ncq_Nrq(2)];
                    %Namefig = ['Fit results of ' flag_row_or_col ' '  num2str(indexq) ' '  CCN2S_struct.TITLEstruct.TITLEstr2 ];
                    Namefig = ['Fit at constant nu = ' num2str(CCN2S_struct.scancq(indexq).scanrq(1).nu,'%10.3e') ' (1/Angstroms)  '  CCN2S_struct.TITLEstruct.TITLEstr1(size(CCN2S_struct.TITLEstruct.TITLEstr1,1),:)];
                case 'col'
                    qtolook = [ 1:CCN2S_struct.Ncq_Nrq(1)];
                    % Namefig = ['Fit results of ' flag_row_or_col ' '  num2str(indexq) ' '  CCN2S_struct.TITLEstruct.TITLEstr2 ];
                    Namefig = ['Fit at constant del = ' num2str(CCN2S_struct.scancq(1).scanrq(indexq).del,'%10.3e') ' (1/Angstroms)  '   CCN2S_struct.TITLEstruct.TITLEstr1(size(CCN2S_struct.TITLEstruct.TITLEstr1,1),:) ];
                    
                otherwise
                    disp('please select rows or cols')
                    return;
            end
            
            
            for iq = qtolook
                
                switch flag_row_or_col
                    
                    case 'row'
                        qvector(iq) = CCN2S_struct.scancq(indexq).scanrq(iq).del;
                        CCNdtV_fit  = CCN2S_struct.scancq(indexq).scanrq(iq).CCNdtV_fit;
                        title_tag = 'col';
                    case 'col'
                        qvector(iq) = CCN2S_struct.scancq(iq).scanrq(indexq).nu;
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
                        XLabelstr = 'qvector\_del in 1/A';
                        
                        switch pp
                            case 3
                                Axislim_vect =  [-5e-4 5e-4 0 2.5e3];%[-60 60 0 5e8];%
                            otherwise
                                Axislim_vect =  [-5e-4 5e-4  -abs(min(fitres_CCN2S.pout(pp).val))-1e-4 abs(max(fitres_CCN2S.pout(pp).val))+1e-4];%[-60 60 0 5e8];%
                        end
                        
                        
                    case 'col'
                        XLabelstr = 'qvector\_nu in 1/A';
                        switch pp
                            case 3
                                Axislim_vect = [-20e-3 20e-3 0 5.0e3];%[-60 60 0 1];% [-60 60 0 5e8];%[-2e-2 2e-2 0 1.5e3];
                            otherwise
                                Axislim_vect = [-20e-3 20e-3 0 1e1];%[-60 60 0 1];%[-60 60 0 5e8];%[-2e-2 2e-2  -abs(min(fitres_CCN2S.pout(pp).val))-1e-4 abs(max(fitres_CCN2S.pout(pp).val))+1e-4];
                        end
                end
                
                DisplayFunctions_XPCS.display_style(fighandle,'figure',Titlestr,Namestr,XLabelstr,YLabelstr,Shading_mode,Axislim_vect,flagPrettyPlot,Colorvector,Colorbarflag,Axisimagestr )
                
            end
            
            %set(gcf,'Name',Namefig);
            
            %
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
                            qvector(iq) = Allscans(iT).CCN2S_struct.scancq(indexq).scanrq(iq).del;
                            pout  = Allscans(iT).CCN2S_struct.scancq(indexq).scanrq(iq).CCNdtV_fit.pout(pp_index);
                            sigp  = Allscans(iT).CCN2S_struct.scancq(indexq).scanrq(iq).CCNdtV_fit.sigp(pp_index);
                            ptitle = Allscans(iT).CCN2S_struct.scancq(indexq).scanrq(iq).CCNdtV_fit.plegend(pp_index).ptitle;
                            Axislim_vect = [1e-5 1e-3 1 1e4];%[1 64 1e-3 1];%[1 64 1 5e8];%
                            XLabelstr = 'qvector\_del in 1/A';
                        case 'col'
                            qvector(iq) = Allscans(iT).CCN2S_struct.scancq(iq).scanrq(indexq).nu;
                            pout  = Allscans(iT).CCN2S_struct.scancq(iq).scanrq(indexq).CCNdtV_fit.pout(pp_index);
                            sigp  = Allscans(iT).CCN2S_struct.scancq(iq).scanrq(indexq).CCNdtV_fit.sigp(pp_index);
                            ptitle = Allscans(iT).CCN2S_struct.scancq(iq).scanrq(indexq).CCNdtV_fit.plegend(pp_index).ptitle;
                            Axislim_vect = [3e-5 3e-2 1 1e4];%[1 64 1e-3 1];%[1 64 1 5e8];%
                            XLabelstr = 'qvector\_nu in 1/A';
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
        
        function [maps2D] = display_2Dmap(CCN2S_struct,fignum,QvalFlag,sim_flag)
            
            [~,~,~,~,~,XCOLlabel,YROWlabel,...
                AXISdet,~,INFOstr] = XPCS_initialize_parameters.TTplot_parameters();
            
            
            [ImageJ,~,SINGLE] = XPCS_initialize_parameters.TTsput_read_ini();
            
            
            Nc = numel(CCN2S_struct.scancq);
            Nr = numel(CCN2S_struct.scancq(1).scanrq);
            
            for cc = 1:Nc
                for rr = 1:Nr
                    maps2D.tau(rr,cc) = CCN2S_struct.scancq(cc).scanrq(rr).CCNdtV_fit.pout(3);
                    maps2D.tau_sigma(rr,cc) = CCN2S_struct.scancq(cc).scanrq(rr).CCNdtV_fit.sigp(3);
                    
                    maps2D.contrast(rr,cc) = CCN2S_struct.scancq(cc).scanrq(rr).CCNdtV_fit.pout(2);
                    maps2D.contrast_sigma(rr,cc) = CCN2S_struct.scancq(cc).scanrq(rr).CCNdtV_fit.sigp(2);
                end
                
            end
            
            %{
            figure(fignum);
            imagesc(maps2D.tau);
            colorbar;
            title('tau');
            caxis([0 1e5]);
            
            figure(fignum+1);
            pcolor(maps2D.tau_sigma);
            colorbar;
            title('sigma of tau');
            caxis([0 1e5]);
            %}
            IIstruct.TITLEstuct = CCN2S_struct.TITLEstruct;
            mask_tau = maps2D.tau>0;
            DisplayTTM_plot.display_II(IIstruct,log10(maps2D.tau.*mask_tau),'Time constants (s) in log scale',QvalFlag,fignum,ImageJ,SINGLE,XCOLlabel,YROWlabel,AXISdet,INFOstr,sim_flag,[1 8]);
            
            DisplayTTM_plot.display_II(IIstruct,maps2D.contrast.*mask_tau,'Contrast',QvalFlag,fignum+1,ImageJ,SINGLE,XCOLlabel,YROWlabel,AXISdet,INFOstr,sim_flag,[1 8]);

            
        end
        
        
        
        %%%%% other functions
        
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
        
        
        
        
    end
end

