classdef FittingFunctions
    % This library contains all the functions which allow us to rotate the
    % detector and the sample
    properties(Constant)
    end
    
    
    methods(Static)
        
       function [fitresult] = fit_2time_corr(CCfunc,fitfunc_paramstr,param_legend,fitfunc_str,fit_range,p_lower,p_upper,p_start)
            % this function fits the time correlation function vs the time
            % and extratcs the time constant. It uses the MATLAB routine
            % fit based on a least square fit.
            %
            % Inputs:
            %   CCfunc:structure containing the different correlation
            %   functions
            %   time_vect: structure containing the different time axis
         
            
            fit_line = fittype(fitfunc_str);
            fitresult_func = str2func([fitfunc_paramstr fitfunc_str]);
            
            ftype=fittype(fit_line); % change this if another function is desired instead

            %Choose the parameter bounds, and starting points for the peak:
            opts=fitoptions(fit_line); % change this if another function is used

            opts.Lower = p_lower; % in alphabetical order [a1 b1]
            opts.Upper = p_upper; 
            opts.StartPoint = p_start; % these are the initial guesses

            for jj=1:numel(CCfunc)
                
                [gfit,gof,fitres] = fit(CCfunc(jj).time_1D(fit_range)',CCfunc(jj).CCNdtV(fit_range)',ftype,opts);
                ci=confint(gfit,.95); % Change .95 value if other that 95% confidence intervals desired
                
               
                fitresult(jj).fitfunc = fitresult_func(CCfunc(jj).time_1D,gfit.a1,gfit.b1,gfit.c1);
                fitresult(jj).x = CCfunc(jj).time_1D(fit_range);
                fitresult(jj).y = CCfunc(jj).CCNdtV(fit_range);
                fitresult(jj).pout = [gfit.a1 gfit.b1 gfit.c1]; 
                fitresult(jj).ci = ci;
                fitresult(jj).plegend = param_legend;
                
            end
            
            
            
        end
        
       function [fitresult] = fit_2time_corr_with_leasqr(CCfunc,fitfunc_str,fit_range,pin,dp, w)
            % this function fits the time correlation function vs the time
            % and extratcs the time constant.
            %
            % Inputs:
            %   CCfunc:structure containing the different correlation
            %   functions
            %   time_vect: structure containing the different time axis
            
            
               
            % fit using leasqrs
            global verbose;
            verbose = [0 1];
            
            opts_struct.stol = 0.000001;
            opts_struct.niter = 100;  
            opts_struct.w = w; % weights (error bar for instance)
            opts_struct.pin=pin; % initial parameters
            opts_struct.dp=dp;  % Fraction to vary for deriv
            
          
            
            for jj=1:numel(CCfunc)
                
                [ycalc,pout,plegend,kvg,iter,corp,covp,covr,stdresid,Z,r2]=leasqr(CCfunc(jj).time_1D(fit_range)',CCfunc(jj).CCNdtV(fit_range),opts_struct.pin,fitfunc_str,opts_struct.stol,opts_struct.niter,sqrt(opts_struct.w),opts_struct.dp);                

                fitresult(jj).fitfunc = ycalc;
                fitresult(jj).x = CCfunc(jj).time_1D(fit_range);
                fitresult(jj).y = CCfunc(jj).CCNdtV(fit_range);
                fitresult(jj).plegend = plegend;
                fitresult(jj).pout = pout;
                sigp = zeros(numel(opts_struct.dp),1);
                sigp(find(opts_struct.dp~=0)) = sqrt(diag(covp));
                fitresult(jj).sigp = sigp;
                fitresult(jj).kvg = kvg;
                fitresult(jj).iter = iter;
                fitresult(jj).corp = corp;
                fitresult(jj).covp = covp;
                fitresult(jj).covr = covr;
                fitresult(jj).stdresid = stdresid;
                fitresult(jj).Z = Z;
                fitresult(jj).r2 = r2;
                
            end
            
            
            
       end
                     
       function [y,legendstruct] = CCN2single_fit(x,p)
           % function TTM_fn(x,p)
           %   for fit of CCNdtV
           %   p(1) = C constant (positive)
           %   p(2) = A_avg
           %   p(3) = tau_avg
           %   p(4) = slope
          
           
           y = p(1) + p(2).*exp(-x./p(3)) +p(4).*x;
           legendstruct(1).ptitle = 'Back.';
           legendstruct(2).ptitle = 'Contrast';
           legendstruct(3).ptitle = 'tau (sec)';
           legendstruct(4).ptitle = 'slope (1/sec)';
           
           
       end
       
       function [fitresult] = fit_tau_with_leasqr(pout_struct,fitfunc_str,fit_range,pin,dp, w)
           % this function fits the time correlation function vs the time
           % and extratcs the time constant.
           %
           % Inputs:
           %   CCfunc:structure containing the different correlation
           %   functions
           %   time_vect: structure containing the different time axis
           
           
           
           % fit using leasqrs
           global verbose;
           verbose = [0 1];
           
           opts_struct.stol = 0.000001;
           opts_struct.niter = 100;
           opts_struct.w = w(fit_range); % weights (error bar for instance)
           opts_struct.pin=pin; % initial parameters
           opts_struct.dp=dp;  % Fraction to vary for deriv
           
           
           [ycalc,pout,plegend,kvg,iter,corp,covp,covr,stdresid,Z,r2]=leasqr(pout_struct.qvector(fit_range)',pout_struct.tau(fit_range)',opts_struct.pin,fitfunc_str,opts_struct.stol,opts_struct.niter,sqrt(opts_struct.w),opts_struct.dp);
           
           fitresult.fitfunc = ycalc;
           fitresult.x = pout_struct.qvector(fit_range);
           fitresult.y = pout_struct.tau(fit_range);
           fitresult.plegend = plegend;
           fitresult.pout = pout;
           fitresult.sigp = sqrt(diag(covp));
           fitresult.kvg = kvg;
           fitresult.iter = iter;
           fitresult.corp = corp;
           fitresult.covp = covp;
           fitresult.covr = covr;
           fitresult.stdresid = stdresid;
           fitresult.Z = Z;
           fitresult.r2 = r2;
           
           
           
           
       end
        
     function [y,legendstruct] = powerlaw_tau(x,p)
           % function TTM_fn(x,p)
           %   for fit of CCNdtV
           %   p(1) = C constant (positive)
           %   p(2) = A_avg
           %   p(3) = tau_avg
          
           
           y = p(1)./x;
           legendstruct(1).ptitle = 'Contrast';
           
           
           
     end
          
     function IInormbb_ref = fit_IInormb_intime(IIbin_struct, N_degree)
           % This function fits the temperal dependence of each pixel in
           % IInormb to a polynomial function whose degree is specified by
           % N_degree
           
            IInormbb = IIbin_struct.IInormbb;
            Xamountb = IIbin_struct.Xamountb;
            
            Nrs = size(IInormbb,1);
            Ncs = size(IInormbb,2);
            Ntb = size(IInormbb,3);
            
            IInormbb_ref = zeros(Nrs,Ncs,Ntb);
            
            for irs = 1:Nrs
                for ics =1:Ncs
                   [p] = polyfit(Xamountb',squeeze(IInormbb(irs,ics,:)),N_degree);
                   [IInormbb_ref(irs,ics,:) ]= polyval(p,Xamountb);
                end
            end
            
           
     end
     
     function IInormbb_ref_avg = fit_IInormb_avg_intime(IInormbb_ref_avg, N_degree)
           % This function fits the temperal dependence of each pixel in
           % IInormb to a polynomial function whose degree is specified by
           % N_degree
           
            %IInormbb = IInormbb_ref_avg.IInormbb;
            Xamountb = IInormbb_ref_avg.scancq(1).scanrq(1).timex;
            
            Nrs = numel(IInormbb_ref_avg.scancq(1).scanrq);
            Ncs = numel(IInormbb_ref_avg.scancq);
            Ntb = numel(Xamountb);
            
            IInormbb_avg_fit = zeros(Nrs,Ncs,Ntb);
            
            for irs = 1:Nrs
                for ics =1:Ncs
                   [p] = polyfit(Xamountb',IInormbb_ref_avg.scancq(ics).scanrq(irs).IInormb_avg,N_degree);
                   [IInormbb_avg_fit(irs,ics,:) ]= polyval(p,Xamountb);
                end
            end
            
            IInormbb_ref_avg.IInormbb_avg_fit= IInormbb_avg_fit;
           
       end 
     
    end
end