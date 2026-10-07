%% BYOM function call_deri.m (calculates the model output)
%
%  Syntax: [Xout,TE,Xout2,zvd] = call_deri(t,par,X0v,glo)
%
% This function calls the explicit function(s) in <simplefun.html
% simplefun.m>. The complete calculation takes place there. As input, it
% gets:
%
% * _t_   the time vector
% * _par_ the parameter structure
% * _X0v_   a vector with initial states and one concentration (scenario number)
% * _glo_  the structure with various types of information (used to be global)
%
% The output _Xout_ provides a matrix with time in rows, and states in
% columns. This function calls <simplefun.html simplefun.m>. The optional
% output _TE_ is not used in this package. The _zvd_ for zero-variate data
% is also not used here. 
%
% This function is for the mixture GUTS model. The external function
% simplefun provides all the state variables.
%
% Author: Tjalling Jager
% Date: December 2021
% Web support: <http://www.debtox.info/byom.html>

%  Copyright (c) 2012-2021, Tjalling Jager.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE-MIT.txt file in the root directory of BYOM. 

function [Xout,TE,Xout2,zvd] = call_deri(t,par,X0v,glo)

% These outputs need to be defined, even if they are not used
Xout2    = []; % additional uni-variate output, not used in this example
zvd      = []; % additional zero-variate output, not used in this example

%% Initial settings
% This package only uses the analytical solution in <simplefun.html
% simplefun.m>. All options for the ODE solver are removed.

% Unpack the vector X0v, which is X0mat for one scenario
X0 = X0v(2:end); % these are the intitial states for a scenario

%% Calculations
% This part calls the explicit model in <simplefun.html simplefun.m> to
% calculate the output (the value of the state variables over time). There
% is generally no need to modify this part. 

c = X0v(1); % the concentration (or scenario number)
t = t(:);   % force t to be a row vector (needed when useode=0)

TE = +inf; % return infinity for time of events
% use an explicit function provided in simplefun!
Xout = simplefun(t,X0,par,c,glo);

%% Output mapping
% _Xout_ contains a row for each state variable. It can be mapped to the
% data. If you need to transform the model values to match the data, do it
% here. 
%
% This set of files is geared towards the use of the analytical solution
% for the damage level and survival in simplefun.

% Note: all calculations take place in simplefun!
