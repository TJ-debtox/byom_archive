%% BYOM, byom_calanus2016.m, reduced GUTS
%
% * Author: Tjalling Jager 
% * Date: December 2021
% * Web support: <http://www.debtox.info/byom.html>
%
% BYOM is a General framework for simulating model systems in terms of
% ordinary differential equations (ODEs). The model itself needs to be
% specified in <derivatives.html derivatives.m> or <simplefun.html
% simplefun.m>, and <call_deri.html call_deri.m> may need to be modified to
% the particular problem as well. The files in the engine directory are
% needed for fitting and plotting. Results are shown on screen but also
% saved to a log file (|results.out|). The files in this directory use an
% analytical solution only, and therefore <derivatives.html derivatives.m>
% will be missing.
%
% *The model:* fitting survival data with the
% <http://www.debtox.info/about_guts.html GUTS> special cases based on the
% (scaled) internal concentration: SD, IT and mixed (or GUTS proper). Here,
% the analytical solution in simplefun.m is used, which is much faster than
% the using the ODEs in derivatives. 
%
% *This script:* C2-naphthalene in _Calanus finmarchicus_. 
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

% observed number of survivors, time in days, conc in uM
DATA{1} = [-1	0	1.03	1.97	4.83	10.3	29.3
    0	56	56	28	28	28	28
    1	56	56	28	28	25	0
    2	56	56	28	28	10	0
    3	56	56	28	26	1	0
    4	56	55	26	21	0	0
    5	56	55	26	18	0	0
    6	56	55	25	6	0	0];
     
%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios.

X0mat(1,:) = [DATA{1}(1,2:end)]; % scenarios (concentrations)
X0mat(2,:) = 1; % initial survival probability
X0mat(3,:) = 0; % initial internal concentration

%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. 

% syntax: par.name = [startvalue fit(0/1) minval maxval optional:log/normal scale (0/1)];

% globals used for GUTS purposes
glo.sel   = 1; % select death mechanism: 1) SD 2) IT
glo.locS  = 1; % location of survival probability in the state variable list
glo.locD  = 2; % location of the dose metric in the state variable list

% syntax: par.name = [startvalue fit(0/1) minval maxval optional:log/normal scale (0/1)];
par.kd = [1  1 1e-3  10  1];   % dominant rate constant, d-1
par.mw = [1  1 0     1e6 1];   % median threshold for survival (uM)
par.hb = [1  1 1e-4  1   1];   % background hazard rate (d-1)
par.bw = [1  1 1e-4  1e6 1];   % killing rate (uM-1 d-1) (SD only)
par.Fs = [1  1 1     100 1];   % fraction spread of threshold distribution (IT only)

switch glo.sel % make sure that right parameters are fitted
    case 1 % for SD ...
        par.Fs(2) = 0; % never fit the threshold spread!
    case 2 % for IT ...
        par.bw(2) = 0; % never fit the killing rate!
end

%% Time vector and labels for plots
% Specify what to calculate and what to plot.

% specify the y-axis labels for each state variable
glo.ylab{1} = 'survival probability';
% specify the x-axis label (same for all states)
glo.xlab    = 'time (days)';
glo.leglab1 = ''; % legend label before the 'scenario' number
glo.leglab2 = [char(181),'M']; % legend label after the 'scenario' number

prelim_checks % script to perform some preliminary checks and set things up
% Note: prelim_checks also fills all the options (opt_...) with defauls, so
% modify options after this call, if needed.

%% Calculations and plotting
% Here, the script is called that will do the calculation and the plotting.
% Options for the plotting can be set using opt_plot (see prelim_checks)

opt_optim.it     = 1; % show iterations of the simplex optimisation (1, default) or not (0)
opt_plot.bw      = 1; % plot in black and white
opt_plot.cn      = 0; % if set to 1, connect model line to points (only for bw=1)
opt_plot.annot   = 2; % annotations in multiplot for fits: 1) box with parameter estimates 2) single legend
opt_plot.statsup = [2]; % vector with states to suppress in plotting fits

par = start_vals_guts(par); % experimental start-value finder; use at your own risk
% Note: start_vals will now overwrite the parameter structure par!

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