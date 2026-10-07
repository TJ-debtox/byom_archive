%% BYOM function simplefun_2.m (the model as explicit equations)
%
%  Syntax: Xout = simplefun_2(t,X0,par,c,glo)
%
% This function calculates the output of the full GUTS model system. Note
% that this function calculates the damage _Di_ analytically. Note that the
% survival probability due to chemical stress is all calculated in
% call_deri.m. This version can deal with time-varying exposure
% concentrations, as long as the exposure profile is specified using
% make_scen as type 2, 3 or 4 (smooth splining requires the ODE version of
% the model in derivatives.m). The analytical solutions are pretty ugly,
% but they seem to work well. However, it does not hurt to compare your
% results to the results of the ODE solver (e.g., by not fitting parameters
% and comparing with glo.use_ode 0 and 1), especially when you do not
% trust the results.
% 
% This function is slightly modified from the standard full model, in that
% _Kiw_ is forced to 1 to produce a double-scaled model. Further, diff_rel
% is set to 1e-3, which was needed to avoid the calculations getting stuck
% in some cases (this is a good idea for the standard GUTS package as well).
%
% NOTE: watch out when the rate constant (_ke_, _kr_, _kc_) get close to each
% other. In the analytical solutions, there are divisions by the difference
% between two rate constants, so they are not allowed to be equal. And,
% there will be numerical problems when they get too close. I try to catch
% these cases, and modify one or more of the rate constants a tiny bit.
% However, it is good to be aware of this limitation (when in doubt: also
% run the ODE version).
% 
% As input, it gets:
%
% * _t_   is the time vector
% * _X0_  is a vector with the initial values for states
% * _par_ is the parameter structure
% * _c_   is the external concentration (or scenario number)
% * _glo_ is the structure with information (normally global)
%
% Time _t_ is handed over as a vector, and scenario name _c_ as single
% number, by <call_deri.html call_deri.m> (you do not have to use them in
% this function). Output _Xout_ (as matrix) provides the output for each
% state at each _t_.
%
% * Author: Tjalling Jager 
% * Date: November 2023
% * Web support: <http://www.debtox.info/byom.html>

% Copyright (c) 2012-2023, Tjalling Jager.
% This source code is licensed under the MIT-style license found in the
% LICENSE-MIT.txt file in the root directory of BYOM. 

%% Start

function Xout = simplefun_2(t,X0,par,c,glo)

%% Unpack initial states
% The state variables enter this function in the vector _X_0. 

% S0  = X0(glo.locS); % survival probability at t=0
Ci0 = X0(glo.locC); % internal concentration at t=0
Di0 = X0(glo.locD); % scaled damage (referenced to internal concentrations) at t=0

%% Unpack parameters
% The parameters enter this function in the structure _par_. The names in
% the structure are the same as those defined in the byom script file.
% The 1 between parentheses is needed as each parameter has 5 associated
% values.

ke   = par.ke(1);   % elimination rate constant
kr   = par.kr(1);   % damage repair rate constant
Kiw  = 1;  % bioconcentration factor (set to 1, so I don't need to change the ugly equations below!)
% mw   = par.mi(1);   % median of threshold distribution (used in call_deri)
% bw   = par.bi(1);   % killing rate (used in call_deri)
% Fs   = par.Fs(1);   % fraction spread of threshold distribution, (-) (used in call_deri)
hb   = par.hb(1);   % background hazard rate

%% Calculate the model output
% This is the actual model, specified as explicit function(s):

S  = exp(-hb*t); % fill S with background mortality for survival
% The calculation of the actual survival probability is done in call_deri
% and only background survival is calculated here. 

Tev = [0 c]; % events setting: without anything else, assume it is constant 
kc  = 0; % assume no disappearance of the text chemical
if glo.timevar == 1 % if we are warned that we have a time-varying concentration ...
    [Tev,kc] = read_scen(-2,c,t,glo); % use read_scen to derive actual exposure concentration
    % the -2 lets read_scen know we are calling from simplefun (and need events and kc)
end

if sum(Tev(:,2)) == 0 && Di0 == 0 && Ci0 == 0 % no need to calculate anything if there is only zero concentrations in Tev
    Ci = zeros(length(t),1);
    Di = zeros(length(t),1);
    Xout(:,[glo.locS glo.locC glo.locD]) = [S Ci Di]; % combine them into a matrix
    return % return as we are done here
end

% initialise internal concentrations and damage with NaNs
Ci = nan(length(t),1);
Di = nan(length(t),1);
Ci(1) = Ci0; % the first element is the starting concentration
Di(1) = Di0; % the first element is the starting concentration
    
diff_rel = 1e-3; % minimum relative difference between the rate constants

if size(Tev,2) == 2 && kc == 0 % then we have a series of constant exposures
    % Even though this is a special case of the static-renewal with
    % degradation below, I still keep it separate. Main reason is that when
    % kc is zero (no degradation) and kr is zero (no repair) this leads to
    % problems (and one of the rate constants needs to be modified). The
    % code below does not have that problem. Still requires ke and kr to be
    % different, though.
    
    if abs(1-ke/kr) < diff_rel % the analytical solution does not allow the two rate constants
        % to be exactly the same (or too close) ...
        ke = ke * (1+diff_rel); % increase ke tiny bit (it should never be zero)
    end
    
    for i = 1:size(Tev,1) % run through all event periods
        
        a = ke * Kiw * Tev(i,2); % modify the uptake fluxe to the new concentration in this period
        % compound parameter a is used in the solution below
        
        if i<size(Tev,1)
            ind_t = (t>Tev(i,1) & t<=Tev(i+1,1)); % find the logical indices for the period
        else
            ind_t = (t>Tev(i,1)); % find the logical indices for the last period
            % if the scenario is longer than the data, this is all zeros,
            % which leads to an empty te, but this does not produce an
            % error so it is fine.
        end
        te = t(ind_t) - Tev(i,1); % find the new part of the time vector, and make sure it starts at zero again
        
        % Thanks to Bob Kooi for using Maple to obtain the solutions
        Ci(ind_t) = (a + exp(-ke * te) * (Ci0 * ke - a)) / ke;
        % use the ugly solution
        Di(ind_t) = -(-ke^2 + ke*kr) / ke / (ke-kr)^2 * (Ci0*kr + Di0*ke - Di0*kr - a) * exp(-kr*te) ...
            - exp(-ke*te) * kr / (ke-kr) * (Ci0*ke - a) / ke - (-a*ke + a*kr) / ke / (ke-kr);

        
        % find new starting values at exact moment of new period
        if i < size(Tev,1) % don't do this for last period
            te = Tev(i+1,1) - Tev(i,1); % start time for new period
            Ci_0_n = (a + exp(-ke * te) * (Ci0 * ke - a)) / ke;
            
            Di0 = -(-ke^2 + ke*kr) / ke / (ke-kr)^2 * (Ci0*kr + Di0*ke - Di0*kr - a) ...
                * exp(-kr*te) - exp(-ke*te) * kr / (ke-kr) * (Ci0*ke - a) / ke ...
                - (-a*ke + a*kr) / ke / (ke-kr);
            
            Ci0 = Ci_0_n; % this is needed as Di0 is also based on old Ci0!
        end
    end
    
elseif size(Tev,2) == 2 && kc > 0 % then we have a static renewal scenario

    % kc , ke and kr cannot be too close to each other The code below
    % attempts to solve that in a clever way (as we have 3 parameters to
    % deal with) by slightly enlarging the smallest difference between two
    % parameters.
    
    k_all  = [kc ke kr]; % place all three rate constant in a vector a
    diff_a = diff([k_all k_all(1)])./k_all; % these are the relative differences between all three
    if all(diff_a<3*diff_rel) % then all three parameters are close together ...
        k_all = k_all - [1 0 -1] .* (diff_rel*k_all); % decrease kc and increase kr a tiny bit
    else
        ind_diff    = abs(diff_a) < diff_rel; % find the difference that is too small
        sign_diff   = sign(diff_a); % and the direction of this difference (so we can enlarge it)
        sign_diff(sign_diff==0) = 1; % if two values are exactly the same, sign will be zero, so make it 1
        k_all(ind_diff) = k_all(ind_diff) - sign_diff(ind_diff) .* (diff_rel*k_all(ind_diff));
    end
    k_all = max(0,k_all); % make sure they don't become negative
    kc    = k_all(1);
    ke    = k_all(2);
    kr    = k_all(3);
    if any(diff([k_all k_all(1)])./k_all) < diff_rel
        error('Several of the rate constants are still too close together ... check what is going on (notify Tjalling ;-) )')
    end
    
    for i = 1:size(Tev,1) % run through all event periods
        
        if i<size(Tev,1)
            ind_t = (t>Tev(i,1) & t<=Tev(i+1,1)); % find the logical indices for the period
        else
            ind_t = (t>Tev(i,1)); % find the logical indices for the last period
            % if the scenario is longer than the data, this is all zeros,
            % which leads to an empty te, but this does not produce an
            % error so it is fine.
        end
        te = t(ind_t) - Tev(i,1); % find the new part of the time vector, and make sure it starts at zero again
        
        a = Tev(i,2); % the new start concentration in this period
        % parameter a is used in the solution below
        
        % Thanks to Bob Kooi for using Maple to obtain the solutions
        Ci(ind_t) = (-0.1e1 / (kc - ke) * Kiw * a * ke * exp(-te * (kc - ke)) ...
            + (Kiw * a * ke + Ci0 * kc - Ci0 * ke) / (kc - ke)) ...
            .* exp(-ke * te);

        % use the ugly solution
        Di(ind_t) = (-kr / (kc - ke) * (Kiw * a * ke / (-kc + kr) * exp(-kc * te + kr * te) ...
            - (Kiw * a * ke + Ci0 * kc - Ci0 * ke) / (kc - ke) * kc / (-ke + kr) ...
            * exp(-ke * te + kr * te) + (Kiw * a * ke + Ci0 * kc - Ci0 * ke) ...
            / (kc - ke) * ke / (-ke + kr) * exp(-ke * te + kr * te)) ...
            + (Kiw * a * ke * kr + kc * Ci0 * kr + kc * Di0 * ke - kc ...
            * Di0 * kr - Ci0 * kr ^ 2 - Di0 * ke * kr + Di0 * kr ^ 2)  ...
            / (ke * kc - kc * kr - ke * kr + kr ^ 2)) .* exp(-kr * te);

        
        % find new starting values at exact moment of new period
        if i < size(Tev,1) % don't do this for last period
            te = Tev(i+1,1) - Tev(i,1); % start time for new period
            Ci_0_n = (-0.1e1 / (kc - ke) * Kiw * a * ke * exp(-te * (kc - ke)) ...
                + (Kiw * a * ke + Ci0 * kc - Ci0 * ke) / (kc - ke)) ...
                * exp(-ke * te);
            
            Di0 = (-kr / (kc - ke) * (Kiw * a * ke / (-kc + kr) * exp(-kc * te + kr * te) ...
                - (Kiw * a * ke + Ci0 * kc - Ci0 * ke) / (kc - ke) * kc / (-ke + kr) ...
                * exp(-ke * te + kr * te) + (Kiw * a * ke + Ci0 * kc - Ci0 * ke) ...
                / (kc - ke) * ke / (-ke + kr) * exp(-ke * te + kr * te)) ...
                + (Kiw * a * ke * kr + kc * Ci0 * kr + kc * Di0 * ke - kc ...
                * Di0 * kr - Ci0 * kr ^ 2 - Di0 * ke * kr + Di0 * kr ^ 2)  ...
                / (ke * kc - kc * kr - ke * kr + kr ^ 2)) * exp(-kr * te);

            Ci0 = Ci_0_n; % this is needed as Di0 is also based on old Ci0!
        end
    end
    
else % we are using linear interpolation in a forcing series
       
    if t(end) > Tev(end,1) % do we ask for more points than in Tev?
        Tev = [Tev; t(end) 0 0]; % add a dummy time point in Tev (from t)
    end

    if abs(1-ke/kr) < diff_rel % the analytical solution does not allow the two rate constants
        % to be exactly the same (or too close) ...
        ke = ke * (1+diff_rel); % increase ke tiny bit
    end
    
    for i = 1:size(Tev,1)-1 % run through all events (not last one that we added)
        
        ind_t = (t>Tev(i,1) & t<=Tev(i+1,1)); % find the logical indices for the period
        te = t(ind_t) - Tev(i,1); % find the new part of the time vector, and make it start at zero again
        
        % Thanks to Bob Kooi for using Maple to obtain the solutions
        a = Tev(i,2); % take as initial concentration the new one in this period
        b = Tev(i,3); % take as slope the new one in this period
        
        Ci(ind_t) = Kiw * b * te + Kiw * a - Kiw / ke * b - exp(-te * ke) ...
            * (Kiw * a * ke - Kiw * b - Ci0 * ke) / ke;
        
        Di(ind_t) = (Kiw * b * ke ^ 2 * te * kr - kr ^ 2 * Kiw ...
            * b * te * ke + Kiw * a * ke ^ 2 * kr - kr ^ 2 * Kiw * a ...
            * ke + kr ^ 2 * (Kiw * a * ke - Kiw * b - Ci0 * ke) ...
            * exp(-kr * te - te * (ke - kr)) - exp(-kr * te) ...
            * ke ^ 2 * (Kiw * a * ke * kr - Kiw * b * ke - Ci0 ...
            * kr ^ 2 - Di0 * ke * kr + Di0 * kr ^ 2) / (ke - kr) ...
            + kr * exp(-kr * te) * (Kiw * a * ke * kr - Kiw * b * ke ...
            - Ci0 * kr ^ 2 - Di0 * ke * kr + Di0 * kr ^ 2) ...
            / (ke - kr) * ke - Kiw * b * ke ^ 2 + Kiw * b * kr ^ 2) ...
            / (ke - kr) / kr / ke;
                
        % find new starting values at exact moment of new period
        if i < size(Tev,1)-1 % don't do this for last period
            te = Tev(i+1,1) - Tev(i,1); % start time for new period
            Ci_0_n = Kiw * b * te + Kiw * a - Kiw / ke * b - exp(-te * ke) ...
                * (Kiw * a * ke - Kiw * b - Ci0 * ke) / ke;
            
            Di0 = (Kiw * b * ke ^ 2 * te * kr - kr ^ 2 * Kiw ...
                * b * te * ke + Kiw * a * ke ^ 2 * kr - kr ^ 2 * Kiw * a ...
                * ke + kr ^ 2 * (Kiw * a * ke - Kiw * b - Ci0 * ke) ...
                * exp(-kr * te - te * (ke - kr)) - exp(-kr * te) ...
                * ke ^ 2 * (Kiw * a * ke * kr - Kiw * b * ke - Ci0 ...
                * kr ^ 2 - Di0 * ke * kr + Di0 * kr ^ 2) / (ke - kr) ...
                + kr * exp(-kr * te) * (Kiw * a * ke * kr - Kiw * b * ke ...
                - Ci0 * kr ^ 2 - Di0 * ke * kr + Di0 * kr ^ 2) ...
                / (ke - kr) * ke - Kiw * b * ke ^ 2 + Kiw * b * kr ^ 2) ...
                / (ke - kr) / kr / ke;

            Ci0 = Ci_0_n; % this is needed as Di0 is also based on old Ci0!
        end
    end
       
end
   
Xout(:,[glo.locS glo.locC glo.locD]) = [S Ci Di]; % combine them into a matrix
