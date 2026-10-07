%% BYOM function simplefun.m (the model as explicit equations)
%
%  Syntax: Xout = simplefun(t,X0,par,c,glo)
%
% This function calculates the output of the model system. It is linked to
% the script files <byom_bioconc_extra.html byom_bioconc_extra.m> and
% <byom_bioconc_start.html byom_bioconc_start.m>. As input, it gets:
%
% * _t_   is the time vector
% * _X0_  is a vector with the initial values for states
% * _par_ is the parameter structure
% * _c_   is the external concentration (or scenario number)
% * _glo_ is the structure with information (normally global)
% 
% Note: _glo_ is now an input to this function rather than a global. This
% makes the code faster (though only very little when using simplefun).
%
% Time _t_ is handed over as a vector, and scenario name _c_ as single
% number, by <call_deri.html call_deri.m> (you do not have to use them in
% this function). Output _Xout_ (as matrix) provides the output for each
% state at each _t_.
%
% * Author: Tjalling Jager 
% * Date: November 2021
% * Web support: <http://www.debtox.info/byom.html>
% * Back to index <walkthrough_byom.html>
% 
%  Copyright (c) 2012-2021, Tjalling Jager.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE-MIT.txt file in the root directory of BYOM. 

%% Start
% Note that, in this example, variable _c_ (the scenario identifier) is not
% used in this function. The treatments differ in their initial exposure
% concentration that is set in X0mat. Also, the structure _glo_ is not used
% here.

function Xout = simplefun(t,X0,par,c,glo)

%% Unpack initial states
% The state variables enter this function in the vector _X_0.

% No initial values used (assumed to be zero)

%% Unpack parameters
% The parameters enter this function in the structure _par_. The names in
% the structure are the same as those defined in the byom script file.
% The 1 between parentheses is needed as each parameter has 5 associated
% values.

ku   = par.ku(1);   % uptake rate constant, L/kg/d-1
ke   = par.ke(1);   % elimination rate constant, d-1

%% Calculate the model output
% This is the actual model, specified as explicit function(s):

ind1 = t<=7;
ind2 = t>7;
Ci = nan(length(t),1);

Ci(ind1) = (ku/ke)*c*(1-exp(-ke*t(ind1))); % internal concentration, accumulation
Ci0      = (ku/ke)*c*(1-exp(-ke*7));       % new initial internal concentration
Ci(ind2) = Ci0 * exp(-ke*(t(ind2)-7));     % internal concentration, depuration

Xout = [Ci]; % combine them into a matrix