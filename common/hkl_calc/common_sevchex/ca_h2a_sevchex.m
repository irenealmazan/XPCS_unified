function [ANGLES,DOCU_anglenames,ANGLES2,DOCU_angle2names,ERROR] = ca_h2a_sevchex(hkl,UB,Energy,params)
% [ANGLES,DOCU_anglenames,ANGLES2,DOCU_angle2names,ERROR] 
%		= calc_h2a_sevchex(hkl,UB,Energy,params)
% Carol Thompson, Brian Stephenson 2014 Oct 25   (modes 0-2,4-7 direct, mode 3 iterative)
%	sevchex assumes CW rotation about z axis (inward)
%	2015-0514 CT edited some parameter names to be consistent with other codes parameter names
%%		params.Alpha and params.Beta field names changed to params.alphatarget and params.betatarget
%%		to be consistent with the ca_h2a of the other diffractometers
%%		added output DOCU_angle2names (which shifted position of ERROR in output)
%
% seven circles: Rho, Nu, Delta, Mu, Eta, Chi, Phi 
% seven heavens: Moon, Mercury, Venus, Sun, Mars, Jupiter, Saturn 
% seven dwarfs: Doc, Grumpy, Happy, Sleepy, Bashful, Sneezy, Dopey
% seven deadly sins:  wrath, greed, sloth, pride, lust, envy, gluttony
% seven virtues: prudence, justice, temperance, courage, faith, hope, charity
% seven hills of Rome: Aventine, Caelian, Capitoline, Esquiline, Palatine, Quirinal, Viminal
% seven sisters colleges: Barnard, Bryn Mawr, Mount Holyoke, Radcliffe, Smith, Vassar, Wellesley 
% seven seals: white horse, red horse, black horse, pale horse, souls of martyrs, earthquake, angels
% seven samarai
% seven seas
% seven wonders of world
% seven horcruxes
% seven year itch
%
% OUTPUT 
%   ANGLES = [Rho, Nu, Delta, Mu, Eta, Chi, Phi] in degrees, one row per hkl
%   DOCU_anglenames = documentation of angle names
%   ANGLES2 = [Alpha, Beta, Two_Theta, PolCorr] degrees, polarization correction (0-1), when 0; exptl intensity zero)
%   DOCU_angle2names = documentation of angles2 names 
%   ERROR = structure of error indicators
%       ERROR.Alpha_min = lowest possible Alpha
%       ERROR.Alpha_max = highest possible Alpha
%       ERROR.Alpha_low = true if Alpha < Alpha_min
%       ERROR.Alpha_hi = true if Alpha > Alpha_max
%   For mode 3:
%       ERROR.BB_lo
%       ERROR.ZZ_hi
%       ERROR.BB2_lo
%       ERROR.ZZ2_hi
%       ERROR.dx = final deviation (from zero) of search minimum, 
%       ERROR.iter = iterations 
%   For modes 4, 5, 6:
%       ERROR.BB3_lo
%       ERROR.ZZ3_hi
%       ERROR.HH1_lo
%
% INPUT 
%   hkl = [h k l] in r.l.u, one row per hkl
%	UB = structure from calc_UB, UB.UB is orientation matrix
%   Energy in eV
%   params is structure of parameters necessary for spectrometer, in this case,
%	params.mode
%	params.sigma, params.tau, params.alphatarget, params.betatarget
%	params.Rho, params.Mu, params.Chi, params.Phi
%
%
%
%   params.mode: integer from 0 to 7
%   all modes have Rho, Mu fixed
%   mode for:        n*x=0   fixed chi
%   ----------------------------------
%   fixed alpha         0       4
%   fixed beta          1       5
%   alpha=beta          2       6
%   fixed phi           3       7
%
%   params.sign_nu controls choice of Hlabx root   (params.sign_nu = +1 or -1)
%   giving sign of nu solution, modes 0-2 and 4-6; usually +
%
%   the following are all in degrees, not all are used in every mode
%   params.sigma and params.tau:	specify surface orientation 
%   params.alphatarget, params.betatarget: specify the target alpha (mode 0,4) or beta (mode 1,5)
%   params.Rho, params.Mu:		fixed values (used in all modes)
%   params.Phi:				fixed values (used for modes 3 and 7)
%   params.Chi:				fixed value (used for mores 4-7)

% Order of output angles for sevchex (best practice is to use same order as from 'pa' command)
% Label with two or more spaces between each name
%    DOCU_anglenames = '[Rho  Nu  Delta  Mu  Eta  Chi  Phi] : sevchex'; 
    DOCU_anglenames = 'Rho  Nu  Delta  Mu  Eta  Chi  Phi'; 
    DOCU_angle2names =  'Alpha  Beta  Two_Theta  PolarizationCorrection';

% for calculations hkl should be in columns
    HKL = hkl';
    UBHKL = UB.UB * [HKL] ./ (2*pi*Energy./fhc); % Matrix of Hphi values, one per column
		
    Rho     = params.Rho;
    Mu      = params.Mu;

    if params.mode == 3 | params.mode == 7
        Phif = params.Phi;
    end

    if params.mode == 4 | params.mode == 5 | params.mode == 6 | params.mode == 7
        Chif = params.Chi;
    end

    Sigma 	= params.sigma;
    Tau		= params.tau;
    nphi = [sind(Sigma).*cosd(Tau); -sind(Sigma).*sind(Tau); cosd(Sigma)];
    
    sign_nu = sign(params.sign_nu);
    
	LEN = length(hkl(:,1));
	ANGLES = NaN*ones(LEN,7);
	ANGLES2 = NaN*ones(LEN,3);
	for ii=1:LEN 
		
		Hphi = UBHKL(:,ii);
        
        normH = norm(Hphi); % Put in test: if zero, exit
		
        C0 = normH^2/2;

        Two_Theta = 2*asind(normH/2);
        
        if params.mode == 0 | params.mode == 1 | params.mode == 2 | params.mode == 4 | params.mode == 5 | params.mode == 6
            
            % For modes 0, 1, 4, 5 get limits for obtainable Alpha:
            
            if params.mode == 0 | params.mode == 1 | params.mode == 4 | params.mode == 5
                Xi = acosd(dot(Hphi,nphi)./normH); % Angle between nphi and Hphi
        
                Alpha_min = Two_Theta/2 - Xi;
                Alpha_max = Two_Theta/2 + Xi;      
                ERROR.Alpha_min(ii) = Alpha_min;
                ERROR.Alpha_max(ii) = Alpha_max;
        
                if params.mode == 0 | params.mode == 4 % Fixed alpha
                    Alpha = params.alphatarget;
                    Beta = asind(dot(Hphi,nphi) - sind(Alpha));
                    ERROR.Alpha_low(ii) = Alpha < Alpha_min;
                    ERROR.Alpha_hi(ii) = Alpha > Alpha_max;
                else % mode == 1, 5 Fixed beta
                    Beta = params.betatarget;
                    Alpha = asind(dot(Hphi,nphi) - sind(Beta));
                    ERROR.Alpha_low(ii) = Alpha < Alpha_min;
                    ERROR.Alpha_hi(ii) = Alpha > Alpha_max;
                end
            else % mode == 2, 6 Alpha = Beta
                Alpha = asind(dot(Hphi,nphi)/2);
                Beta = Alpha;
                ERROR.Alpha_low(ii) = 0;
                ERROR.Alpha_hi(ii) = 0;
            end
            
            if ERROR.Alpha_low(ii) | ERROR.Alpha_hi(ii)
            
                Nu = NaN;
                Delta = NaN;
                Eta = NaN;
                Chi = NaN;
                Phi = NaN;
                
            else
        
                % Calculate Eta, Chi, Phi rotation matrix from known
                % transformations of H and n between the mu and phi frames
                % using Busing and Levy procedure
                
                Hlab3 = dot(Hphi,nphi)./cosd(Alpha) - C0.*tand(Alpha);
                Hlab = [sign_nu.*sqrt(2*C0 - C0.^2 - Hlab3.^2); ... % positive root gives pos Nu
                        -C0; ...
                        Hlab3];
            
                Hphi_norm = Hphi./normH; % This will fail at Q=0
            
                Tphi3 = cross(Hphi_norm,nphi);
                Tphi3 = Tphi3./norm(Tphi3);
            
                Tphi2 = cross(Tphi3,Hphi_norm);
                Tphi2 = Tphi2./norm(Tphi2);
            
                Tphi = [Hphi_norm Tphi2 Tphi3];
           
                % Get Hmu, nmu
            
                nlab = [0; -sind(Alpha); cosd(Alpha)]; 
            
                MR = [cosd(Mu)       -sind(Mu).*cosd(Rho)  -sind(Mu).*sind(Rho); ...
                      sind(Mu)        cosd(Mu).*cosd(Rho)   cosd(Mu).*sind(Rho); ...
                      0                        -sind(Rho)             cosd(Rho)];
            
                Hmu = MR*Hlab;
                nmu = MR*nlab;
            
                Hmu_norm = Hmu./normH;
            
                Tmu3 = cross(Hmu_norm,nmu);
                Tmu3 = Tmu3./norm(Tmu3);
            
                Tmu2 = cross(Tmu3,Hmu_norm);
                Tmu2 = Tmu2./norm(Tmu2);
            
                Tmu = [Hmu_norm Tmu2 Tmu3];
            
                R = Tphi*(Tmu');
            
                % This matrix satisfies Hphi = R*Hmu, nphi = R*nmu
            
                % Solve for Eta, Chi, Phi from R = R_phi*R_chi*R_eta (rot matrices)
                % Below works for cosd(Chi) not near zero
            
                Chi = asind(R(3,1));
                Eta = atan2d(-R(3,2),R(3,3));
                Phi = atan2d(R(2,1),R(1,1));
            
                if params.mode == 0 | params.mode == 1 | params.mode == 2     
                
                    % calc of Nu and Delta, then done:
                
                    Delta = asind((Hlab3 + tand(Rho).*(C0 - 1))./(tand(Rho).*sind(Rho) + cosd(Rho)));
        
                    sin_nu = Hlab(1)./cosd(Delta);
                    cos_nu = (1 - C0 + sind(Delta).*sind(Rho))./(cosd(Delta).*cosd(Rho));
                    Nu = atan2d(sin_nu, cos_nu); % May need error check for not real inputs
                
                else % params.mode = 4, 5, 6
                    
                    % Rotate solution around kin, 
                    % Sy is kin in Phi frame
                    % must be the same for both sets of Eta, Chi, Phi
                    % find new Eta, Phi that match Chif
                    
                    Sy = [(sind(Eta).*sind(Phi) + sind(Chi).*cosd(Eta).*cosd(Phi)).*sind(Rho) ...
                        + ((sind(Chi).*sind(Eta).*cosd(Mu) - cosd(Chi).*sind(Mu)).*cosd(Phi) ...
                        - cosd(Eta).*cosd(Mu).*sind(Phi)).*cosd(Rho); ...
                        (-sind(Eta).*cosd(Phi) + sind(Chi).*cosd(Eta).*sind(Phi)).*sind(Rho) ...
                        + ((sind(Chi).*sind(Eta).*cosd(Mu) - cosd(Chi).*sind(Mu)).*sind(Phi) ...
                        + cosd(Eta).*cosd(Mu).*cosd(Phi)).*cosd(Rho); ...
                          (-cosd(Chi).*sind(Eta).*cosd(Mu) - sind(Chi).*sind(Mu) ).*cosd(Rho) ...
                        - cosd(Chi).*cosd(Eta).*sind(Rho)];
                    
                    Chi = Chif;
                    
                    % Use Sy(3) to solve for Eta
                    
                    % AA*cos(x) + BB*sin(x) = CC = A*sin(x+y)
                    % y = atan(AA/BB)
                    % x = asin(CC*cos(y)/BB) - y
                    
                    AA3 = -cosd(Chif).*sind(Rho);
                    BB3 = -cosd(Chif).*cosd(Mu).*cosd(Rho);
                    CC3 = Sy(3) + sind(Chif).*sind(Mu).*cosd(Rho);
                   
                    if (abs(BB3) < 1.e-20) % e.g. Rho = 0 & Mu = 90 & Chif = 0
                        ERROR.BB3_lo(ii) = 1;
                        Eta = NaN;
                        Phi = NaN;
                    else
                        ERROR.BB3_lo(ii) = 0;
                        YY3 = atand(AA3./BB3); % Keep YY3 near zero
                        ZZ3 = CC3.*cosd(YY3)./BB3;
                        if abs(ZZ3) > 1
                            ERROR.ZZ3_hi(ii) = 1;
                            Eta = NaN;
                            Phi = NaN;
                        else
                            ERROR.ZZ3_hi(ii) = 0;
                            Eta = asind(ZZ3) - YY3; % Keep Eta near zero
                        end
                            
                        % Now solve for Phi
                    
                        AA1 = sind(Chif).*cosd(Eta).*sind(Rho) + (sind(Chif).*sind(Eta).*cosd(Mu) - cosd(Chif).*sind(Mu)).*cosd(Rho);
                        BB1 = sind(Eta).*sind(Rho) - cosd(Eta).*cosd(Mu).*cosd(Rho);
                           
                        if (hypot(AA1,BB1) < 1.e-20)
                            ERROR.HH1_lo(ii) = 1;
                            Phi = NaN;
                        else
                            ERROR.HH1_lo(ii) = 0;
                            Phi = atan2d(BB1.*Sy(1) + AA1.*Sy(2), AA1.*Sy(1) - BB1.*Sy(2));
                        end
                      
                    end
                   
                end        
                
            end
            
        elseif params.mode == 3 % dot(n,x) = 0, fixed Phi
            
            % Iterative solution for Eta and Chi, converges immediately for
            % Mu = 0, 180
            % First calc nchi, Hchi in chi frame
            
            Phi = Phif;
            
            P_inv = [cosd(Phi) sind(Phi) 0; ...
                    -sind(Phi) cosd(Phi) 0; ...
                        0       0        1];
            Hchi = P_inv * Hphi;
            nchi = P_inv * nphi;
            
            % Below works for Mu near zero or 180; 
            % could swap which dot products are used for chi and phi 
            % for Mu near 90 or 270
            
            Eta0 = 0; % Initial guess for Eta
            
            for jj = [1:100] % Iterate until Eta converges
                
                % First solve for Chi to get dot(n,x) = 0
                % R is R_chi * R_eta * R_mu * R_rho (rot matrices)
                
                %Rx = [sind(Chi).*sind(Eta).*sind(Mu) + cosd(Chi).*cosd(Mu); ...
                %      cosd(Eta).*sind(Mu); ...
                %     -cosd(Chi).*sind(Eta).*sind(Mu) + sind(Chi).*cosd(Mu)];
                %
                % Split into coefficients of cosd(Chi), sind(Chi), constant
            
                Rxc = [cosd(Mu); 0; -sind(Eta0).*sind(Mu)];
                Rxs = [sind(Eta0).*sind(Mu); 0; cosd(Mu)];
                Rx0 = [0; cosd(Eta0).*sind(Mu); 0];
                
                AA = dot(nchi,Rxc);
                BB = dot(nchi,Rxs);
                CC = -dot(nchi,Rx0);
            
                % AA*cos(x) + BB*sin(x) = CC = A*sin(x+y)
                % y = atan(AA/BB)
                % x = asin(CC*cos(y)/BB) - y
                   
                if (abs(BB) < 1.e-20) % e.g. Sigma = 0 & Mu = 90
                    ERROR.BB_lo(ii) = 1;
                    Eta = NaN;
                    Chi = NaN;
                    break;
                else
                    ERROR.BB_lo(ii) = 0;
                    YY = atand(AA./BB); % Keep YY near zero
                    ZZ = CC.*cosd(YY)./BB;
                    if abs(ZZ) > 1
                        ERROR.ZZ_hi(ii) = 1;
                        Eta = NaN;
                        Chi = NaN;
                        break;
                    else
                        ERROR.ZZ_hi(ii) = 0;
                        Chi = asind(ZZ) - YY; % Keep Chi near zero
                    end
                end
                
                % Next solve for Eta to get dot(H,y) = -C0
                
                %Ry = [sind(Chi).*cosd(Eta).*sind(Rho) ...
                %   + (sind(Chi).*sind(Eta).*cosd(Mu) - cosd(Chi).*sind(Mu)).*cosd(Rho); ...
                %    -sind(Eta).*sind(Rho) + cosd(Eta).*cosd(Mu).*cosd(Rho); ...
                %    (-cosd(Chi).*sind(Eta).*cosd(Mu) - sind(Chi).*sind(Mu)).*cosd(Rho) ...
                %   - cosd(Chi).*cosd(Eta).*sind(Rho)];
                %
                % Split into coefficients of cosd(Eta), sind(Eta), constant
                
                Ryc = [sind(Chi).*sind(Rho); cosd(Mu).*cosd(Rho); -cosd(Chi).*sind(Rho)];
                Rys = [sind(Chi).*cosd(Mu).*cosd(Rho); -sind(Rho); -cosd(Chi).*cosd(Mu).*cosd(Rho)];
                Ry0 = [-cosd(Chi).*sind(Mu).*cosd(Rho); 0; -sind(Chi).*sind(Mu).*cosd(Rho)];
                
                AA2 = dot(Hchi,Ryc);
                BB2 = dot(Hchi,Rys);
                CC2 = -C0 -dot(Hchi,Ry0);
            
                % AA*cos(x) + BB*sin(x) = CC = A*sin(x+y)
                % y = atan(AA/BB)
                % x = asin(CC*cos(y)/BB) - y
                   
                if (abs(BB2) < 1.e-20) % e.g. Rho = 0 & Mu = 90
                    ERROR.BB2_lo(ii) = 1;
                    Eta = NaN;
                    Chi = NaN;
                    break;
                else
                    ERROR.BB2_lo(ii) = 0;
                    YY2 = atand(AA2./BB2); % Keep YY near zero
                    ZZ2 = CC2.*cosd(YY2)./BB2;
                    if abs(ZZ2) > 1
                        ERROR.ZZ2_hi(ii) = 1;
                        Eta = NaN;
                        Chi = NaN;
                        break;
                    else
                        ERROR.ZZ2_hi(ii) = 0;
                        Eta = asind(ZZ2) - YY2; % Keep Eta near zero
                    end
                end
               
                % Check for convergence
                if (abs(Eta - Eta0)) <= 1.e-3
                    break
                else
                    Eta0 = Eta;
                end     
            end
            ERROR.iter(ii) = jj;
            ERROR.dx(ii) = Eta - Eta0;
            
        else % params.mode == 7 Fixed Chi, Phi
            
            Chi = Chif;
            Phi = Phif;
            
            % Split Sy into coefficients of cosd(Eta), sind(Eta), and constant
            
            Syc = [sind(Chif).*cosd(Phif).*sind(Rho) - cosd(Mu).*sind(Phif).*cosd(Rho); ...
                   sind(Chif).*sind(Phif).*sind(Rho) + cosd(Mu).*cosd(Phif).*cosd(Rho); ...
                  -cosd(Chif).*sind(Rho)];
            
            Sys = [sind(Phif).*sind(Rho) + sind(Chif).*cosd(Mu).*cosd(Phif).*cosd(Rho); ...
                  -cosd(Phif).*sind(Rho) + sind(Chif).*cosd(Mu).*sind(Phif).*cosd(Rho); ...
                  -cosd(Chif).*cosd(Mu).*cosd(Rho)];
            
            Sy0 = [-cosd(Chif).*sind(Mu).*cosd(Phif).*cosd(Rho); ...
                   -cosd(Chif).*sind(Mu).*sind(Phif).*cosd(Rho); ...
                   -sind(Chif).*sind(Mu).*cosd(Rho)];
            
            AA = dot(Hphi,Syc);
            BB = dot(Hphi,Sys);
            CC = -C0 - dot(Hphi,Sy0);
            
            % AA*cos(x) + BB*sin(x) = CC = A*sin(x+y)
            % y = atan(AA/BB)
            % x = asin(CC*cos(y)/BB) - y
                   
            if (abs(BB) < 1.e-20) % e.g. Rho = 0 & Mu = 90 & Chif = 0
                ERROR.BB_lo(ii) = 1;
                Eta = NaN;
            else
                ERROR.BB_lo(ii) = 0;
                YY = atand(AA./BB); % Keep YY near zero
                Eta = asind(CC.*cosd(YY)./BB) - YY; % Keep Eta near zero, need error check here
            end
        
        end

        % For modes 3 to 7, calculate Delta, Nu
        
        if params.mode > 2
            % Columns of S matrix:
            Sx = [(sind(Chi).*sind(Eta).*sind(Mu) + cosd(Chi).*cosd(Mu)).*cosd(Phi) ...
                - cosd(Eta).*sind(Mu).*sind(Phi); ...
                  (sind(Chi).*sind(Eta).*sind(Mu) + cosd(Chi).*cosd(Mu)).*sind(Phi) ...
                + cosd(Eta).*sind(Mu).*cosd(Phi); ...
                  -cosd(Chi).*sind(Eta).*sind(Mu) + sind(Chi).*cosd(Mu)];    

            Sz = [(-sind(Eta).*sind(Phi) - sind(Chi).*cosd(Eta).*cosd(Phi)).*cosd(Rho) ...
                + ((sind(Chi).*sind(Eta).*cosd(Mu) - cosd(Chi).*sind(Mu)).*cosd(Phi) ...
                - cosd(Eta).*cosd(Mu).*sind(Phi)).*sind(Rho); ...
                   (sind(Eta).*cosd(Phi) - sind(Chi).*cosd(Eta).*sind(Phi)).*cosd(Rho) ...
                + ((sind(Chi).*sind(Eta).*cosd(Mu) - cosd(Chi).*sind(Mu)).*sind(Phi) ...
                + cosd(Eta).*cosd(Mu).*cosd(Phi)).*sind(Rho); ...
                   (-cosd(Chi).*sind(Eta).*cosd(Mu) - sind(Chi).*sind(Mu) ).*sind(Rho) ...
                + cosd(Chi).*cosd(Eta).*cosd(Rho)];

            Hlabx = dot(Sx,Hphi);
            Hlabz = dot(Sz,Hphi);
            
            Delta = asind((Hlabz + tand(Rho).*(C0 - 1))./(tand(Rho).*sind(Rho) + cosd(Rho)));
            sin_nu = Hlabx./cosd(Delta);   
            cos_nu = (1 - C0 + sind(Delta).*sind(Rho))./(cosd(Delta).*cosd(Rho));
            if ~isreal(sin_nu) | ~isreal(cos_nu) % May not need this check
                ERROR.sin_nu(ii) = sin_nu;
                ERROR.cos_nu(ii) = cos_nu;
                Nu = NaN;
            else
                Nu = atan2d(sin_nu, cos_nu); 
            end
        end
        
        % For modes 3 and 7, calculate Alpha, Beta
        
        if params.mode == 3 | params.mode == 7
            % Column of S matrix:
            
            Sy = [(sind(Eta).*sind(Phi) + sind(Chi).*cosd(Eta).*cosd(Phi)).*sind(Rho) ...
                + ((sind(Chi).*sind(Eta).*cosd(Mu) - cosd(Chi).*sind(Mu)).*cosd(Phi) ...
                - cosd(Eta).*cosd(Mu).*sind(Phi)).*cosd(Rho); ...
                 (-sind(Eta).*cosd(Phi) + sind(Chi).*cosd(Eta).*sind(Phi)).*sind(Rho) ...
                + ((sind(Chi).*sind(Eta).*cosd(Mu) - cosd(Chi).*sind(Mu)).*sind(Phi) ...
                + cosd(Eta).*cosd(Mu).*cosd(Phi)).*cosd(Rho); ...
                  (-cosd(Chi).*sind(Eta).*cosd(Mu) - sind(Chi).*sind(Mu) ).*cosd(Rho) ...
                - cosd(Chi).*cosd(Eta).*sind(Rho)];

            C1 = - dot(nphi,Sy);
            Alpha = asind(C1);
            Beta = asind(dot(Hphi,nphi) - C1);
            ERROR.Alpha_low(ii) = 0;
            ERROR.Alpha_hi(ii) = 0;
            
        end

	% (PolarizationCorr*100 is percent (Theoretical*Corr = Experimental Intensity)	
	PolarizationCorr = (cosd(Delta).*cosd(Nu)).^2 + sind(Delta).^2;
	
	ANGLES2(ii,1) = Alpha;
        ANGLES2(ii,2) = Beta;
        ANGLES2(ii,3) = Two_Theta;
	ANGLES2(ii,4) = PolarizationCorr;

        ANGLES(ii,:) = [Rho, Nu, Delta, Mu, Eta, Chi, Phi];
        
	end
	
end

% Helper functions

% Given the angles calculated from these algorithms, and sigma tau, 
% What is alpha and beta  (to either compare with those assumed  

function [Alpha,Beta] = calc_alpha(ANGLES,params)
    Sigma 	= params.sigma;
    Tau		= params.tau;
    nphi = [sind(Sigma).*cosd(Tau); -sind(Sigma).*sind(Tau); cosd(Sigma)];
    
    %['[Rho Nu Delta Mu Eta Chi Phi]'];
    Rho	= ANGLES(:,1);
    Mu 	= ANGLES(:,4);
    Eta = ANGLES(:,5);
    Phi = ANGLES(:,7);
   
    Sy = [(sind(Eta).*sind(Phi) + sind(Chi).*cosd(Eta).*cosd(Phi)).*sind(Rho) ...
          + ((sind(Chi).*sind(Eta).*cosd(Mu) - cosd(Chi).*sind(Mu)).*cosd(Phi) ...
          - cosd(Eta).*cosd(Mu).*sind(Phi)).*cosd(Rho); ...
           (-sind(Eta).*cosd(Phi) + sind(Chi).*cosd(Eta).*sind(Phi)).*sind(Rho) ...
           + ((sind(Chi).*sind(Eta).*cosd(Mu) - cosd(Chi).*sind(Mu)).*sind(Phi) ...
           + cosd(Eta).*cosd(Mu).*cosd(Phi)).*cosd(Rho); ...
             (-cosd(Chi).*sind(Eta).*cosd(Mu) - sind(Chi).*sind(Mu) ).*cosd(Rho) ...
           - cosd(Chi).*cosd(Eta).*sind(Rho)];

    C1 = - dot(nphi,Sy);
    Alpha = asind(C1);
    Beta = asind(dot(Hphi,nphi) - C1);

end



