classdef FittingFunctions
    % This library contains all the functions which allow us to rotate the
    % detector and the sample
    properties(Constant)
    end
    
    
    methods(Static)
        
        
        function IInormbb_ref = fit_IInormb_intime(IIbin_struct, N_degree)
            % This function fits the temporal dependence of each pixel in
            % IInormb to a polynomial function whose degree is specified by
            % N_degree
            
            IInormbb = IIbin_struct.II;
            timeX = IIbin_struct.timeX;
            
            Nrs = size(IInormbb,1);
            Ncs = size(IInormbb,2);
            Ntb = size(IInormbb,3);
            
            IInormbb_ref = zeros(Nrs,Ncs,Ntb);
            
            for irs = 1:Nrs
                for ics =1:Ncs
                    [p] = polyfit(timeX',squeeze(IInormbb(irs,ics,:)),N_degree);
                    [IInormbb_ref(irs,ics,:) ]= polyval(p,timeX);
                end
            end
            
            
        end
        
        function IInormbb_ref_avg = fit_IInormb_avg_intime(IInormbb_ref_avg, N_degree)
            % This function fits the temperal dependence of each pixel in
            % IInormb to a polynomial function whose degree is specified by
            % N_degree
            
            %IInormbb = IInormbb_ref_avg.IInormbb;
            timeXb = IInormbb_ref_avg.scancq(1).scanrq(1).timex;
            
            Nrs = numel(IInormbb_ref_avg.scancq(1).scanrq);
            Ncs = numel(IInormbb_ref_avg.scancq);
            Ntb = numel(timeXb);
            
            IInormbb_avg_fit = zeros(Nrs,Ncs,Ntb);
            
            for irs = 1:Nrs
                for ics =1:Ncs
                    [p] = polyfit(timeXb',IInormbb_ref_avg.scancq(ics).scanrq(irs).IInormb_avg,N_degree);
                    [IInormbb_avg_fit(irs,ics,:) ]= polyval(p,timeXb);
                end
            end
            
            IInormbb_ref_avg.IInormbb_avg_fit= IInormbb_avg_fit;
            
        end        
        
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
            %   p(4) = slope; we force to have a decreasing background
            %   p(5) =
            
            
            %y = p(1) + p(2).*exp(-x./p(3)) +p(4).*x;
            %y = p(1) + p(2).*exp(-x./p(3)) - abs(p(4)).*x;
            y = p(1) + p(2).*exp(-x./p(3)) - p(4)^2.*x;
            
            legendstruct(1).ptitle = 'Back.';
            legendstruct(2).ptitle = 'Contrast';
            legendstruct(3).ptitle = 'tau (sec)';
            legendstruct(4).ptitle = 'slope (1/sec)';
            
            
        end
      
        
          function [y,legendstruct] = CCN2single_fit_double_exp(x,p)
            % function TTM_fn(x,p)
            %   for fit of CCNdtV
            %   p(1) = C constant (positive)
            %   p(2) = A_avg
            %   p(3) = tau_avg
            %   p(4) = slope; we force to have a decreasing background
            %   p(5) =
            
            
            %y = p(1) + p(2).*exp(-x./p(3)) +p(4).*x;
            %y = p(1) + p(2).*exp(-x./p(3)) - abs(p(4)).*x;
            y = p(1) + p(2).*exp(-x./p(3))+p(4).*exp(-x./p(5)) - p(6)^2.*x;
            
            legendstruct(1).ptitle = 'Back.';
            legendstruct(2).ptitle = 'Contrast fast';
            legendstruct(3).ptitle = 'tau fast (sec)';
            legendstruct(4).ptitle = 'Contrast slow';
            legendstruct(5).ptitle = 'tau slow (sec)';
            legendstruct(6).ptitle = 'slope (1/sec)';
            
            
        end
      
        
        function [pin_array] = prepare_coefs_forglobalfit(pin_iTT,iT,qvector)
            
            num_coefs =  size(pin_iTT,2);
            num_Qs = numel(qvector);
            pin_array = zeros(1,num_coefs*num_Qs);
            
            for qq = 1:num_Qs
                
                for ii = [1 3 4]
                    pin_array((qq-1)*num_coefs+ii) = pin_iTT(iT,ii); % constant background
                    % pin_array((qq-1)*num_coefs+3) = pin_iTT(iT,3); % time constant
                    %pin_array((qq-1)*num_coefs+4) = pin_iTT(iT,4); % slope background
                end
                % contrast
                pin_array((qq-1)*num_coefs+2) = pin_iTT(iT,2)./sqrt(1+(2*pi./(1e3*qvector(qq))).^-2);%pin_iTT(iT,2)./sqrt(1+(2*pi./(1e3*qvector(qq))).^-2);
                
            end
            
        end
        
        
        function [fitresult] = fit_2time_corr_with_leasqr_globalfit(CCfunc,fitfunc_str,fit_range,pin,dp, w)
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
            
            x_array = [];
            y_array = [];
            %pin_array = [];
            dp_array = [];
            w_array = [];
            for jj=1:numel(CCfunc)
                x_array = [x_array CCfunc(jj).time_1D(fit_range)];
                y_array = [y_array CCfunc(jj).CCNdtV(fit_range)];
                w_array = [w_array w'];
                %pin_array = [pin_array pin];
                dp_array = [dp_array dp];
            end
            
            num_coefs = round(numel(dp_array)/numel(CCfunc)); % number of coefs per data set
            
            opts_struct.stol = 0.000001;
            opts_struct.niter = 100;
            opts_struct.w = w_array; % weights (error bar for instance)
            opts_struct.pin = pin; % initial parameters (already prepared in an array by prep_coefs_forglobalfit
            opts_struct.dp = dp_array;  % Fraction to vary for deriv
            
            [ycalc,pout,plegend,kvg,iter,corp,covp,covr,stdresid,Z,r2]=leasqr(x_array',y_array,opts_struct.pin,fitfunc_str,opts_struct.stol,opts_struct.niter,sqrt(opts_struct.w),opts_struct.dp);
            
            for jj=1:numel(CCfunc)
                fitresult(jj).fitfunc = ycalc((jj-1)*numel(fit_range)+1:(jj)*numel(fit_range));
                fitresult(jj).x = CCfunc(jj).time_1D(fit_range);
                fitresult(jj).y = CCfunc(jj).CCNdtV(fit_range);
                fitresult(jj).plegend = plegend;
                fitresult(jj).pout = pout((jj-1)*num_coefs+1:(jj)*num_coefs);
                sigp = zeros(numel(opts_struct.dp),1);
                sigp(find(opts_struct.dp~=0)) = sqrt(diag(covp));
                fitresult(jj).sigp = sigp((jj-1)*num_coefs+1:(jj)*num_coefs);
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
               
        
        function [y,legendstruct] = CCN2single_globalfit(x,p)
            % This function does a global fit of the entire data set of one
            % time correlation function. In particular, p(2) follows a
            % function of the momentum transfer. x is time.
            %   for fit of CCNdtV
            %   p(1) = C constant (positive)
            %   p(2) = A_avg
            %   p(3) = tau_avg
            %   p(4) = slope; we force to have a decreasing background
     
            
            num_coefs = 4; % number of coeficients per data set
           % counter_p = 1; % counter running through all the coefficients for all the data stets
            num_data_sets = round((numel(p)/num_coefs)); % number of data sets which is equal to the number of Qs to fit simultaneously
            x_data_set_length = round(numel(x)/num_data_sets);
            
            y_conc = [];
            
            for kk = 1:num_data_sets
                % prepare the x values for a data set
                x_data_set = x((kk-1)*x_data_set_length+1:kk*x_data_set_length);
                p_one_data_set = zeros(1,num_coefs);
                
                % prepare coefficients for a single data set
                for ii = 1:num_coefs
                    index_p = (kk-1)*num_coefs + ii;
                    p_one_data_set(ii) = p(index_p);
                   
                   
                    
                   % counter_p = counter_p +1;
                end
                
                % calculate the fitted value for a data set
                y_one_data_set = p_one_data_set(1)...
                    + p_one_data_set(2).*exp(-x_data_set./p_one_data_set(3)) ...
                    - p_one_data_set(4)^2.*x_data_set;
                 
                % concatenate the guessed values for the output
                y_conc = [y_conc y_one_data_set'];
                
            end
            
            y = y_conc';
            
            legendstruct(1).ptitle = 'Back. Global';
            legendstruct(2).ptitle = 'Contrast Global';
            legendstruct(3).ptitle = 'tau (sec) Global';
            legendstruct(4).ptitle = 'slope (1/sec) Global';
            
            
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
        
        function [y,legendstruct] = invsquarelaw_tau(x,p)
            % function TTM_fn(x,p)
            %   for fit of CCNdtV
            %   p(1) = C constant (positive)
            %   p(2) = A_avg
            %   p(3) = tau_avg
            
            
            y = p(1)./x.^2;
            legendstruct(1).ptitle = 'Offset';
            
            
            
        end
        
        
    end
end