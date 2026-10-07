function dX = deri(t,X,c,par)

% This function is used by script_guts.m to generate test data for GUTS.
% 
% Author     : Tjalling Jager (tjalling@debtox.info)
% Date       : April 2017
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2023, Tjalling Jager.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE-MIT.txt file in the root directory of BYOM. 

kd  = par(1); % elimination (dominant) rate constant
mw  = par(2); % threshold
bw  = par(3); % killing rate
hb  = par(4); % background hazard

Dw = X(1); % scaled damage level
H  = X(2); % cumulative hazard

dDw = kd * (c - Dw); % change in scaled damage
dH  = hb + bw * max(0,(Dw - mw)); % hazard for survival, based on scaled damage

dX = [dDw;dH];