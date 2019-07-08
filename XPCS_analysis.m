
classdef XPCS_analysis
    % This library contains all the functions which allow us to analyze the
    % time correlatin functions
    properties(Constant)
    end
    
    
    methods(Static)
  
        function Qval_struct = calculate_qval(XCEN_II,YCEN_II,Ncq_vect,Nrq_vect,sim_flag)
            
            if sim_flag
                Qval_struct.nu = (Nrq_vect-YCEN_II ); %nu direction, in 1/Angstroms
                Qval_struct.del = (Ncq_vect-XCEN_II);% del direction  in 1/Angstroms
                
            else
                [D_ds,kvector,lambda,pixel_size,th_Bragg,SplitQ_CTR,t_osc] = XPCS_initialize_parameters.TT_experiment();
                Qval_struct.nu = (kvector*(Nrq_vect-YCEN_II )*pixel_size/D_ds)*1e-10; % nu direction, in 1/Angstroms
                Qval_struct.del = (kvector*(Ncq_vect-XCEN_II)*pixel_size/D_ds)*sind(th_Bragg)*1e-10; % del direction  in 1/Angstroms
                
            end
            
        end
        
        
        
        function [itt_range_struct] = prepare_subrange_for2corr(iT,Xamount,XCENV,YCENV,XWIDV,YWIDV,tminv,tmaxv)
            % subrange for 2-time correlation function calculation
            
            itt_range_struct.xx = [-XWIDV(iT):XWIDV(iT)]; % Pixel locations relative to center 
            itt_range_struct.yy = [-YWIDV(iT):YWIDV(iT)]; % Pixel locations relative to center
            itt_range_struct.ixtt = XCENV(iT) + itt_range_struct.xx;
            itt_range_struct.iytt = YCENV(iT) + itt_range_struct.yy;
            itt_range_struct.Ncs = length(itt_range_struct.xx);
            itt_range_struct.Nrs = length(itt_range_struct.yy);
            itt_range_struct.ittt = find(Xamount > tminv(iT) & Xamount < tmaxv(iT));
            itt_range_struct.Nts = length(itt_range_struct.ittt);
            
            
        end
        
      
        function [IIbin_struct] = bin_scans_in_time(Read_Singlescan_struct,Singlescan_struct,iT,tbinV,ImageJ,POINTSUMS)
            
            % Use sub-range of pixels and time
            
            iytt = Singlescan_struct.itt_range_struct.iytt;
            ixtt = Singlescan_struct.itt_range_struct.ixtt;
            ittt = Singlescan_struct.itt_range_struct.ittt;
            
          
            
            IInormbs = Read_Singlescan_struct.IIstruct.IInormb(iytt,ixtt,ittt);
            
            Nrs = size(IInormbs,1);
            Ncs = size(IInormbs,2);
            Nts = size(IInormbs,3);
            
            timestampXs = Read_Singlescan_struct.IIstruct.timestampX(ittt);
            Xamounts = Read_Singlescan_struct.IIstruct.Xamount(ittt);
            Xsteps = Read_Singlescan_struct.IIstruct.Xsteps(ittt);
            
            tbin = tbinV(iT);
            
            % First bin scans in time
            if tbin > 1
                tbin = floor(tbin);
                Ntb = floor(Nts/tbin);
                IIbin_struct.Ntb =  Ntb;
                IInormbe = reshape(IInormbs(:,:,1:Ntb*tbin),Nrs,Ncs,tbin,Ntb);
                IIbin_struct.IInormbb = squeeze(sum(IInormbe,3));
                Xamounte = reshape(Xamounts(1:Ntb*tbin),tbin,Ntb);
                Xstepse = reshape(Xsteps(1:Ntb*tbin),tbin,Ntb);
                IIbin_struct.Xamountb = squeeze(mean(Xamounte,1));
                IIbin_struct.Xstepsb = squeeze(mean(Xstepse,1));
                timestampe = reshape(timestampXs(1:Ntb*tbin),tbin,Ntb);
                IIbin_struct.timeb = squeeze(mean(timestampe,1));
                
                if ~isempty(POINTSUMS)
                    POINTSUMSB(:,1) = floor((POINTSUMS(:,1)-1+ImageJ)/tbin)+1-ImageJ;
                    POINTSUMSB(:,2) = floor((POINTSUMS(:,2)+ImageJ)/tbin) - ImageJ;
                else
                    POINTSUMSB = POINTSUMS;
                end
                
            else
                Ntb = Nts;
                IIbin_struct.Ntb = Nts;
                IIbin_struct.IInormbb = IInormbs;
                IIbin_struct.Xamountb = Xamounts';
                IIbin_struct.Xstepsb = Xsteps;
                IIbin_struct.timeb = timestampXs;
                POINTSUMSB = POINTSUMS;
            end
            
            if isempty(POINTSUMSB)
                IIbin_struct.POINTSB=[1 IIbin_struct.Ntb]-ImageJ;
            else
                IIbin_struct.POINTSB = POINTSUMSB;
            end
                        
            IIbin_struct.YROWpts = [1:Nrs] - ImageJ;
            IIbin_struct.XCOLpts = [1:Ncs] - ImageJ;
            IIbin_struct.TITLEstuct = Read_Singlescan_struct.IIstruct.TITLEstuct;
        end
        
        function [IInormbb_ref_avg] = calc_IInormb_avg_inpixels(IIbin_struct,iT,hwttr_allT,hwttc_allT,wrq_allT,wcq_allT,offsetcc_allT,offsetrc_allT,sim_flag,mask)
           % This function fits the temperal dependence of each pixel in
           % IInormb to a polynomial function whose degree is specified by
           % N_degree
           
            IInormbb = IIbin_struct.IInormbb;
            Xamountb = IIbin_struct.Xamountb;
            
            Nrs = size(IInormbb,1);
            Ncs = size(IInormbb,2);
            Ntb = size(IInormbb,3);
            
            hwttr = hwttr_allT(iT);
            hwttc = hwttc_allT(iT);
            wrq = wrq_allT(iT);
            wcq = wcq_allT(iT);
            offsetcc = offsetcc_allT(iT);
            offsetrc = offsetrc_allT(iT);
            
            mask_struct.offsetr = mask.offsetrc_allT(iT);
            mask_struct.hwttr = ones(1,numel([-hwttr:hwttr]));
            if mask.hwttr_allT(iT) ~= 0
                midindex_hwttr = round(numel(mask_struct.hwttr)/2);%+mask_struct.offsetr;
                mask_struct.hwttr(midindex_hwttr-mask.hwttr_allT(iT):midindex_hwttr+mask.hwttr_allT(iT)) =  zeros(numel([-mask.hwttr_allT(iT):mask.hwttr_allT(iT)]),1);
            end
            
            mask_struct.offsetc = mask.offsetcc_allT(iT);
            mask_struct.hwttc = ones(1,numel([-hwttc:hwttc]));
            if mask.hwttc_allT(iT) ~= 0
                midindex_hwttc = round(numel(mask_struct.hwttc)/2);%+mask_struct.offsetc;
                mask_struct.hwttc(midindex_hwttc-mask.hwttc_allT(iT):midindex_hwttc+mask.hwttc_allT(iT)) = zeros(numel([-mask.hwttc_allT(iT):mask.hwttc_allT(iT)]),1);
            end
            
            
            
            
           
            
            ittccen =1 + (Ncs - 1)/2 ;% index of col center
            ittrcen = 1 + (Nrs - 1)/2; % index of row center        

            
            Ncq = 2*wcq + 1;
            Nrq = 2*wrq + 1;
         
            for icq = 1:Ncq
                 offttc = (icq-wcq-1)*(2*hwttc(icq)+1)+ offsetcc;  
                %indexc_non_masc = [-hwttc:hwttc] - mask_struct.hwttc;
                ittc_index = round(ittccen + offttc + [-hwttc:hwttc]) ;
                ittc = ittc_index(find(mask_struct.hwttc~=0));
               
                for irq = 1:Nrq
                    offttr = (irq - wrq - 1)*(2*hwttr+1)+ offsetrc;
                    ittr = round(ittrcen + offttr + [-hwttr:hwttr]);
                    
                    Qval_struct = XPCS_analysis.calculate_qval(ittccen,ittrcen,ittccen + offttc,ittrcen + offttr,sim_flag);
                    
                    IInormbb_ref_avg.scancq(icq).scanrq(irq).IInormb_avg = squeeze(mean(mean(IInormbb(ittr,ittc,:),1),2));
                    
                    IInormbb_ref_avg.qvector.nu(irq) = Qval_struct.nu;
                    
                    IInormbb_ref_avg.scancq(icq).scanrq(irq).timex = Xamountb;
                end
                IInormbb_ref_avg.qvector.del(icq) = Qval_struct.del;
            end
            
            IInormbb_ref_avg.Ncq_Nrq = [Ncq Nrq];
            IInormbb_ref_avg.TITLEstruct = IIbin_struct.TITLEstuct;

       end
        
        function [CCN2_struct,IIbin_struct] = calc_2time_corr(IIbin_struct,iT,N_degree,flag_mean_or_poly)
                     
            IInormbb =  IIbin_struct.IInormbb;
            
            Nrs = size(IInormbb,1);
            Ncs = size(IInormbb,2);
            Ntb = size(IInormbb,3);
            
            switch flag_mean_or_poly
                case 'mean'
                    IInormbb_ref = mean(IInormbb,3);
                case 'poly'
                    IInormbb_ref = FittingFunctions.fit_IInormb_intime(IIbin_struct, N_degree(iT));
                    IIbin_struct.N_degree = N_degree(iT);
            end
            
            IIbin_struct.IInormbb_ref = IInormbb_ref;
            
            
            dI = IInormbb - IInormbb_ref;
            dlnI = IInormbb./IInormbb_ref - 1;
             
            % Calc 2-time using time average mean, no ensemble
            
            %CCN2_struct.IIM2 = NaN*ones(Nrs,Ncs,Ntb,Ntb);
            CCN2_struct.CCN2 = NaN*ones(Nrs,Ncs,Ntb,Ntb);
            
            for ii = 1:Ntb
                for jj = 1:ii
                    %CCN2_struct.IIM2(:,:,ii,jj) = dI(:,:,ii).*dI(:,:,jj);
                    %CCN2_struct.IIM2(:,:,jj,ii) = CCN2_struct.IIM2(:,:,ii,jj);
                   
                    CCN2_struct.CCN2(:,:,ii,jj) = dlnI(:,:,ii).*dlnI(:,:,jj);
                    CCN2_struct.CCN2(:,:,jj,ii) = CCN2_struct.CCN2(:,:,ii,jj);
                end
            end
            %IID2 = diag(IIM2);
            %CC2 = IIM2./sqrt(IID2*IID2'); % Normalized to make diagonal unity
            CCN2_struct.TITLEstuct = IIbin_struct.TITLEstuct;
        end
        
        function [CCN2avg_struct,qvector_struct] = from_CCN2V_to_CCN2avg(Single_scan_struct,iT,hwttr_allT,hwttc_allT,wrq_allT,wcq_allT,offsetcc_allT,offsetrc_allT,sim_flag)
            
            hwttr = hwttr_allT(iT);
            hwttc = hwttc_allT(iT);
            wrq = wrq_allT(iT);
            wcq = wcq_allT(iT);
            offsetcc = offsetcc_allT(iT);
            offsetrc = offsetrc_allT(iT);
            
            
            Nts = size(Single_scan_struct.CCN2_struct.CCN2,3);
            Ncs = size(Single_scan_struct.CCN2_struct.CCN2,2);
            Nrs = size(Single_scan_struct.CCN2_struct.CCN2,1);
            
            ittccen =1 + (Ncs - 1)/2 ;% index of col center
            ittrcen = 1 + (Nrs - 1)/2; % index of row center        

            CCN2V = Single_scan_struct.CCN2_struct.CCN2;
            
            Ncq = 2*wcq + 1;
            Nrq = 2*wrq + 1;
         
            for icq = 1:Ncq
                offttc = (icq-wcq-1)*(2*hwttc+1)+ offsetcc;               
                ittc = round(ittccen + offttc + [-hwttc:hwttc]) ;
               
                for irq = 1:Nrq
                    offttr = (irq - wrq - 1)*(2*hwttr+1)+ offsetrc;
                    ittr = round(ittrcen + offttr + [-hwttr:hwttr]);
                    
                    Qval_struct = XPCS_analysis.calculate_qval(ittccen,ittrcen,ittccen + offttc,ittrcen + offttr,sim_flag);

                   
                    CCN2avg_struct.scancq(icq).scanrq(irq).CCN2avg = squeeze(mean(mean(CCN2V(ittr,ittc,:,:),1),2));
                    CCN2avg_struct.scancq(icq).scanrq(irq).timex = Single_scan_struct.IIbin_struct.Xamountb;
                    CCN2avg_struct.scancq(icq).scanrq(irq).TITLEstr2V = Single_scan_struct.IIbin_struct.TITLEstuct.TITLEstr2;
                    CCN2avg_struct.scancq(icq).scanrq(irq).nu = Qval_struct.nu;
                    CCN2avg_struct.scancq(icq).scanrq(irq).del = Qval_struct.del;
                    CCN2avg_struct.qvector.nu(irq) = Qval_struct.nu;
                    CCN2avg_struct.boxcenterrc.offttr(irq) = ittrcen + offttr;
                end
               CCN2avg_struct.qvector.del(icq) = Qval_struct.del;
               CCN2avg_struct.boxcenterrc.offttc(icq) = ittccen + offttc;
            end
            CCN2avg_struct.Nts_Ncs_Nrs = [Nts Ncs Nrs];
            CCN2avg_struct.Ncq_Nrq = [Ncq Nrq];
            CCN2avg_struct.TITLEstruct = Single_scan_struct.CCN2_struct.TITLEstuct;
            CCN2avg_struct.ittccen = ittccen ;
            CCN2avg_struct.ittrcen = ittrcen ;
            CCN2avg_struct.hwttr = hwttr;
            CCN2avg_struct.hwttc = hwttc;
            CCN2avg_struct.wrq = wrq;
            CCN2avg_struct.wcq = wcq;
            CCN2avg_struct.offsetcc = offsetcc;
            CCN2avg_struct.offsetrc = offsetrc;
            
            qvector_struct = CCN2avg_struct.qvector;
        end
       
        function [CCN2avg_struct,qvector_struct] = from_CCN2V_to_CCN2avg_several_regions(Single_scan_struct,iT,hwttr_allT,hwttc_allT,wrq_allT,wcq_allT,offsetcc_allT,offsetrc_allT,D_ds,kvector,pixel_size,th_Bragg,mask)
            
            hwttr = hwttr_allT(iT);
            hwttc = hwttc_allT(iT);
            wrq = wrq_allT(iT);
            wcq = wcq_allT(iT);
            offsetcc = offsetcc_allT(iT);
            offsetrc = offsetrc_allT(iT);
            
            mask_struct.offsetr = mask.offsetrc_allT(iT);
            mask_struct.hwttr = ones(1,numel([-hwttr:hwttr]));
            if mask.hwttr_allT(iT) ~= 0
                midindex_hwttr = round(numel(mask_struct.hwttr)/2);%+mask_struct.offsetr;
                mask_struct.hwttr(midindex_hwttr-mask.hwttr_allT(iT):midindex_hwttr+mask.hwttr_allT(iT)) =  zeros(numel([-mask.hwttr_allT(iT):mask.hwttr_allT(iT)]),1);
            end
            
            mask_struct.offsetc = mask.offsetcc_allT(iT);
            mask_struct.hwttc = ones(1,numel([-hwttc:hwttc]));
            if mask.hwttc_allT(iT) ~= 0
                midindex_hwttc = round(numel(mask_struct.hwttc)/2);%+mask_struct.offsetc;
                mask_struct.hwttc(midindex_hwttc-mask.hwttc_allT(iT):midindex_hwttc+mask.hwttc_allT(iT)) = zeros(numel([-mask.hwttc_allT(iT):mask.hwttc_allT(iT)]),1);
            end
            
            
            
            
            
            Nts = size(Single_scan_struct.CCN2_struct.CCN2,3);
            Ncs = size(Single_scan_struct.CCN2_struct.CCN2,2);
            Nrs = size(Single_scan_struct.CCN2_struct.CCN2,1);
            
            ittccen =1 + (Ncs - 1)/2 ;% index of col center
            ittrcen = 1 + (Nrs - 1)/2; % index of row center        

            CCN2V = Single_scan_struct.CCN2_struct.CCN2;
            
            Ncq = 2*wcq + 1;
            Nrq = 2*wrq + 1;
      
            for icq = 1:Ncq
                offttc = (icq-wcq-1)*(2*hwttc(icq)+1)+ offsetcc;  
                ittc_index = round(ittccen + offttc + [-hwttc:hwttc]) ;
                ittc = ittc_index(find(mask_struct.hwttc~=0));
                
                for irq = 1:Nrq
                    offttr = (irq - wrq - 1)*(2*hwttr+1)+ offsetrc;
                    %indexr_non_masc = [-hwttr:hwttr] - mask_struct.hwttr;
                    %ittr_index = round(ittrcen + offttr + [-hwttr:hwttr]) ;
                    %ittr = ittr_index(find( mask_struct.hwttr~=0));
                    ittr = round(ittrcen + offttr + [-hwttr:hwttr]);
                    
                    Qval_struct = XPCS_analysis.calculate_qval(ittccen,ittrcen,ittccen + offttc,ittrcen + offttr,D_ds,kvector,pixel_size,th_Bragg);

                   
                    CCN2avg_struct.scancq(icq).scanrq(irq).CCN2avg = squeeze(mean(mean(CCN2V(ittr,ittc,:,:),1),2));
                    CCN2avg_struct.scancq(icq).scanrq(irq).timex = Single_scan_struct.IIbin_struct.Xamountb;
                    CCN2avg_struct.scancq(icq).scanrq(irq).TITLEstr2V = Single_scan_struct.IIbin_struct.TITLEstuct.TITLEstr2;
                    CCN2avg_struct.scancq(icq).scanrq(irq).nu = Qval_struct.nu;
                    CCN2avg_struct.scancq(icq).scanrq(irq).del = Qval_struct.del;
                    CCN2avg_struct.qvector.nu(irq) = Qval_struct.nu;
                    CCN2avg_struct.boxcenterrc.offttr(irq) = ittrcen + offttr;
                end
               CCN2avg_struct.qvector.del(icq) = Qval_struct.del;
               CCN2avg_struct.boxcenterrc.offttc(icq) = ittccen + offttc;
            end
            CCN2avg_struct.Nts_Ncs_Nrs = [Nts Ncs Nrs];
            CCN2avg_struct.Ncq_Nrq = [Ncq Nrq];
            CCN2avg_struct.TITLEstruct = Single_scan_struct.CCN2_struct.TITLEstuct;
            CCN2avg_struct.ittccen = ittccen ;
            CCN2avg_struct.ittrcen = ittrcen ;
            CCN2avg_struct.hwttr = hwttr;
            CCN2avg_struct.hwttc = hwttc;
            CCN2avg_struct.wrq = wrq;
            CCN2avg_struct.wcq = wcq;
            CCN2avg_struct.offsetcc = offsetcc;
            CCN2avg_struct.offsetrc = offsetrc;
            CCN2avg_struct.mask = mask_struct;
            
            qvector_struct = CCN2avg_struct.qvector;
        end
       
        
        function [CCN2S_struct] = from_CCN2avg_to_CCN2S(CCN2avg_struct)
            % This function builds the CCN_structure where we store the
            % integrated 2-times correlation function. The inputs are:
            %   Nrq: the number of pixels/number of q's
            %   CCN_struct: the structure containing the integrated 2-times correlation function
            %   TITLEstr2V: the title string which specifies the name of
            %   the scans in iiT
            
            
            
            Ndt = CCN2avg_struct.Nts_Ncs_Nrs(1);
            idt = 1:Ndt;
            
            
            Ncq = CCN2avg_struct.Ncq_Nrq(1);
            Nrq = CCN2avg_struct.Ncq_Nrq(2);
            
            for icq = 1:Ncq
                
                CCNdt = NaN*ones(Nrq,Ndt);
                
                for irq = 1:Nrq
                    time_scan = CCN2avg_struct.scancq(icq).scanrq(irq).timex - CCN2avg_struct.scancq(icq).scanrq(irq).timex(1); % Delta in time from region start
                    CCN2S = CCN2avg_struct.scancq(icq).scanrq(irq).CCN2avg(idt,idt);
                    for ii = idt
                        CCNdt(irq,ii) = mean(diag(CCN2S,ii-1));
                    end
                    CCN2S_struct.scancq(icq).scanrq(irq).CCNdtV = CCNdt(irq,:);
                    CCN2S_struct.scancq(icq).scanrq(irq).time_1D = time_scan;
                    CCN2S_struct.scancq(icq).scanrq(irq).nu = CCN2avg_struct.scancq(icq).scanrq(irq).nu;
                    CCN2S_struct.scancq(icq).scanrq(irq).del = CCN2avg_struct.scancq(icq).scanrq(irq).del;
                end
            end
            CCN2S_struct.Ncq_Nrq = CCN2avg_struct.Ncq_Nrq;
            CCN2S_struct.Ndt = Ndt;
            CCN2S_struct.ittccen= CCN2avg_struct.ittccen;
            CCN2S_struct.ittrcen =  CCN2avg_struct.ittccen ;
            
            CCN2S_struct.TITLEstruct = CCN2avg_struct.TITLEstruct;
        end
        
        function [CCN2Sfit_struct] = fit_CCN2S(CCN2S_struct,fitfunc_param_str,fitfunc_str,param_legend,fit_range,fig_plot_fit_results)
            
            % initialize fitting results
            a1 = zeros(numel(CCN2S_struct),numel(CCN2S_struct(1).scanq),1);
            b1 = zeros(numel(CCN2S_struct),numel(CCN2S_struct(1).scanq),1);
            
          
            figh = figure(fig_plot_fit_results);
            clf;
            
            for iT = 1:numel(CCN2S_struct)
                for irq = 1:numel(CCN2S_struct(iT).scanq)
                    p_lower =[0 0 0]; 
                    p_upper =[7E-1 100 1e-2]; 
                    p_start = [6e-2 50 1e-3];
                    
                    [CCN2Sfit_struct(iT).scanq(irq).CCNdtV_fit] = FittingFunctions.fit_2time_corr(CCN2S_struct(iT).scanq(irq),fitfunc_param_str,param_legend,fitfunc_str,fit_range,p_lower,p_upper,p_start);                              

                    %DisplayFunctions_XPCS.display_fit_result(CCN_struct,iT,irq,figh);
                    
                end
               
             CCN2Sfit_struct(iT).TITLEstruct = CCN2S_struct.TITLEstruct;                        
            end
        end
        
        function [CCN2S_struct] = fit_CCN2S_with_leasqr(CCN2S_struct,iT,pin_iiT,dp_iiT,fitfunc_str,fit_range,indexq,flag_row_or_col,fig_plot_fit_results)
            
            figure(fig_plot_fit_results);
            clf;
            
            switch flag_row_or_col
                
                case 'row'
                    qtolook = [ 1:CCN2S_struct.Ncq_Nrq(2)];
                case 'col'
                    qtolook = [ 1:CCN2S_struct.Ncq_Nrq(1)];
                otherwise
                    disp('please select rows or cols')
                    return;
            end
            
            for iq = qtolook
                
                
                switch flag_row_or_col
                    
                    case 'row'
                        CCN2S_struct_singlescan = CCN2S_struct.scancq(indexq).scanrq(iq);
                        
                    case 'col'
                        CCN2S_struct_singlescan = CCN2S_struct.scancq(iq).scanrq(indexq);
                        
                    otherwise
                        disp('please select rows or cols')
                        return;
                end
                
                pin = pin_iiT(iT,:);
                dp =  dp_iiT(iT,:);
                w = ones(length(fit_range),1);
                
                [ fitres] = FittingFunctions.fit_2time_corr_with_leasqr(CCN2S_struct_singlescan,fitfunc_str,fit_range,pin,dp, w);
                
                %                     CCN2Sfit_struct(iT).TITLEstruct = CCN2S_struct(iT).TITLEstruct;
                %
                switch flag_row_or_col
                    case 'row'
                        CCN2S_struct.scancq(indexq).scanrq(iq).CCNdtV_fit =  fitres;
                    case 'col'
                        CCN2S_struct.scancq(iq).scanrq(indexq).CCNdtV_fit =  fitres;
                    otherwise
                        disp('please select rows or cols')
                        return;
                end
                
                
                
            end
            
            %DisplayFunctions_XPCS.display_fit_result(CCN_struct,iT,figh);
            
            
        end
        
        function [pout_struct] = fit_tau_with_leasqr(pout_struct,fitfunc_str,qrange)
                      
            
            for iT = 1:numel(pout_struct)  
                
                    pout_struct_fitrange.qvector = pout_struct(iT).qvector(pout_struct(iT).fitrange);
                    sigma = pout_struct(iT).sigma(pout_struct(iT).fitrange);
                    pout_struct_fitrange.tau = pout_struct(iT).tau(pout_struct(iT).fitrange);
                    [~,low_index] = find(pout_struct_fitrange.qvector>qrange(1));
                    fit_range = low_index;
                
                    pin = [1];
                    dp =  [1]'*0.0001;
                    w = 1./sigma;
                    
                    [ pout_struct(iT).pout] = FittingFunctions.fit_tau_with_leasqr(pout_struct_fitrange,fitfunc_str,fit_range,pin,dp, w);
                
                   
            end
        end
        
        function [pout_avg_struct] = calc_average_row_col(pcolor_struct,ppvector)
            
            for pp = ppvector
             
                pout = pcolor_struct(pp).pout;
                sigma = pcolor_struct(pp).sigma;
                
             
                 
                pout_avg_struct(pp).weight = 1./sigma.^2;
                pout_avg_struct(pp).weight_pout = pout.*pout_avg_struct(pp).weight;
                
                pout_avg_struct(pp).sum_row = sum(pout_avg_struct(pp).weight_pout,1,'omitnan');
                pout_avg_struct(pp).sum_col = sum(pout_avg_struct(pp).weight_pout,2,'omitnan');
                
                pout_avg_struct(pp).sum_row_weight = sum(pout_avg_struct(pp).weight,1,'omitnan');
                pout_avg_struct(pp).sum_col_weight = sum(pout_avg_struct(pp).weight,2,'omitnan');
                
                pout_avg_struct(pp).weight_mean_row = pout_avg_struct(pp).sum_row./pout_avg_struct(pp).sum_row_weight;
                pout_avg_struct(pp).weight_mean_col = pout_avg_struct(pp).sum_col./pout_avg_struct(pp).sum_col_weight;
                
                pout_avg_struct(pp).weight_mean = mean(pout_avg_struct(pp).weight_mean_row);
             
                
                pout_avg_struct(pp).del = pcolor_struct(pp).del;
                pout_avg_struct(pp).nu = pcolor_struct(pp).nu;
                pout_avg_struct(pp).ptitle = pcolor_struct(pp).ptitle;
                
                pout_avg_struct(pp).TITLEstruct = pcolor_struct(pp).TITLEstruct;
            end
           
        end
            
        function [fft_sumROIs,freq_array] = calc_fft_sumROIs(Allscans,Read_Allscans,time_index,ROI_index,fignum)  
            % this function calculates the fft of the integrated intensity
            % in the ROIS and plots the spectrum
            
            time_array = [flipud(-Read_Allscans.IIstruct.timeX(time_index) );Read_Allscans.IIstruct.timeX(time_index) ];
            
            delta_freq = 2*pi/(time_array(end)-time_array(1));
            
            freq_array = [-numel(time_array)*delta_freq/2:delta_freq:-delta_freq 0:delta_freq:(numel(time_array)/2-1)*delta_freq];
            
            fft_sumROIs = fftshift(fft([Allscans.ROIS_struct.sum.sumROIS(time_index,ROI_index);flipud(Allscans.ROIS_struct.sum.sumROIS(time_index,ROI_index))]));
            
            figure(fignum);
            clf
            
            subplot(121);
            plot(Read_Allscans.IIstruct.timeX(time_index),Allscans.ROIS_struct.sum.sumROIS(time_index,ROI_index),'r','LineWidth',3.0);
            xlabel('time (s)');
             set(gca,'FontSize',30);
            
            subplot(122);
            plot(freq_array,fft_sumROIs,'LineWidth',3.0);
            title(['ROIS #' num2str(ROI_index)]);
            xlabel('freq (Hz)');
            set(gca,'FontSize',30);
        end
        
    end
    
end