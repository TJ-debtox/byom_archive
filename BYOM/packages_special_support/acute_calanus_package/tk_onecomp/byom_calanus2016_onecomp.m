%% BYOM, byom_calanus_2016_onecomp.m 
%
% * Author: Tjalling Jager 
% * Date: December 2021
% * Web support: <http://www.debtox.info/byom.html>
%
% BYOM is a General framework for simulating model systems in terms of
% ordinary differential equations (ODEs). The model itself needs to be
% specified in <derivatives.html derivatives.m>, and <call_deri.html
% call_deri.m> may need to be modified to the particular problem as well.
% The files in the engine directory are needed for fitting and plotting.
% Results are shown on screen but also saved to a log file (results.out).
%
% *The model:* An organism is exposed to a chemical in its surrounding
% medium. The animal accumulates the chemical according to standard
% one-compartment first-order kinetics.  
%
% *This script:* One-compartment TK of C2-naphthalene in _Calanus finmarchicus_.
% 
%  Copyright (c) 2012-2021, Tjalling Jager.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE-MIT.txt file in the root directory of BYOM. 

%% Initial things
% Make sure that this script is in a directory somewhere *below* the BYOM
% folder.

clear, clear global % clear the workspace and globals
global DATA W X0mat % make the data set and initial states global variables
global glo          % allow for global parameters in structure glo
diary off           % turn of the diary function (if it is accidentaly on)
% set(0,'DefaultFigureWindowStyle','docked'); % collect all figure into one window with tab controls
set(0,'DefaultFigureWindowStyle','normal'); % separate figure windows

pathdefine(0) % set path to the BYOM/engine directory (option 1 uses parallel toolbox)
glo.basenm  = mfilename; % remember the filename for THIS file for the plots
glo.saveplt = 0; % save all plots as (1) Matlab figures, (2) JPEG file or (3) PDF (see all_options.txt)

%% The data set
% Data are entered in matrix form, time in rows, scenarios (exposure
% concentrations) in columns. First column are the exposure times, first
% row are the concentrations or scenario numbers. The number in the top
% left of the matrix indicates how to calculate the likelihood:
%
% * -1 for multinomial likelihood (for survival data)
% * 0  for log-transform the data, then normal likelihood
% * 0.5 for square-root transform the data, then normal likelihood
% * 1  for no transformation of the data, then normal likelihood

DATA{1} = [0.5	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2
    0	2.634920501	1.128226067	2.067572648	1.239467263	2.133021081	3.968671089	2.232029925	1.71563171
    0.25	346.903685	326.4266006	330.6922294	NaN	NaN	NaN	NaN	NaN
    0.5	616.6576112	592.2330082	677.762585	NaN	NaN	NaN	NaN	NaN
    0.75	762.6419853	750.075518	854.8605579	840.4228865	NaN	NaN	NaN	NaN
    1	1187.964104	1224.36273	1075.951767	1121.508449	NaN	NaN	NaN	NaN
    2	2027.432656	2037.744655	1638.995127	NaN	NaN	NaN	NaN	NaN
    3	2524.564775	2471.508177	2675.329348	NaN	NaN	NaN	NaN	NaN
    4	2685.458367	2782.061779	2345.599143	3347.038305	NaN	NaN	NaN	NaN
    4.25	2194.05948	2119.389788	2501.947601	3182.207062	NaN	NaN	NaN	NaN
    4.5	2309.995772	2353.777302	2290.190177	2226.329017	NaN	NaN	NaN	NaN
    4.75	2462.489785	2316.686876	2141.947838	1862.544099	NaN	NaN	NaN	NaN
    5	1772.724311	2170.229973	2019.317227	1786.46631	NaN	NaN	NaN	NaN
    6	2076.314829	1390.663864	742.5848544	NaN	NaN	NaN	NaN	NaN
    7	1800.811532	1862.208204	1223.001179	1050.560266	NaN	NaN	NaN	NaN
    8	1920.780515	1490.778938	1228.435888	1602.794556	NaN	NaN	NaN	NaN];

W{1} = 21 * ones(size(DATA{1})-1); % each point is a pooled sample of 21 animals

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios.

X0mat = [1.2    % the scenarios (here nominal concentrations) 
         0];  % initial values state 1 (internal concentrations)

%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. 

% syntax: par.name = [startvalue fit(0/1) minval maxval];
par.ke    = [0.2    1 0.01 100 1];  % elimination rate constant, d-1
par.ku    = [900    1 0.01 1e6 1];  % uptake rate constant, L/kg/d

%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% used, based on the data set

% specify the y-axis labels for each state variable
glo.ylab{1} = ['internal concentration (',char(181),'mol/kg)'];
% specify the x-axis label (same for all states)
glo.xlab    = 'time (days)';
glo.leglab1 = ''; % legend label before the 'scenario' number
glo.leglab2 = [char(181),'M']; % legend label after the 'scenario' number

prelim_checks % script to perform some preliminary checks and set things up
% Note: prelim_checks also fills all the options (opt_...) with defauls, so
% modify options after this call, if needed.

%% Calculations and plotting
% Here, the function is called that will do the calculation and the plotting.
% Options for the plotting can be set using opt_plot (see prelim_checks.m).
% Options for the optimsation routine can be set using opt_optim. Options
% for the ODE solver are part of the global glo. 

opt_optim.it = 1; % show iterations of the simplex optimisation (1, default) or not (0)
opt_plot.bw  = 1; % plot in black and white
opt_plot.cn  = 0; % if set to 1, connect model line to points (only for bw=1)
glo.useode   = 0; % use the analytical solution in simplefun.m

% optimise and plot (fitted parameters in par_out)
par_out = calc_optim(par,opt_optim); % start the optimisation
calc_and_plot(par_out,opt_plot); % calculate model lines and plot them

%% Profiling the likelihood
% By profiling you make robust confidence intervals for one or more of your
% parameters. Use the names of the parameters as they occurs in your
% parameter structure _par_ above. This can be a single string (e.g.,
% 'kd'), a cell array of strings (e.g., {'kd','ke'}), or 'all' to profile
% all fitted parameters. 
% 
% Options for profiling can be set using opt_prof (see prelim_checks.m).

opt_prof.detail   = 2; % detailed (1) or a coarse (2) calculation
opt_prof.subopt   = 10; % number of sub-optimisations to perform to increase robustness

% UNCOMMENT LINE(S) TO CALCULATE
par_better = calc_proflik(par_out,{'all'},opt_prof,opt_optim);  % calculate a profile
if ~isempty(par_better)                 % if the profiling found a better optimum ...
    calc_and_plot(par_better,opt_plot); % calculate model lines and plot them
end
