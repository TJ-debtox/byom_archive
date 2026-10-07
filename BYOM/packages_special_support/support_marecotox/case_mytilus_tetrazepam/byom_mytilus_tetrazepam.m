%% BYOM, byom_mytilus_tetrazepam.m
% This is the file to reproduce the case study in Section 3.3.3 of the
% book (Fig. 3.3).
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
% one-compartment first-order kinetics, specified by an uptake rate
% constant (_ku_) and an elimination rate (_ke_). 
%
% *This script:* Data from Gomez et al. (2012),
% <http://dx.doi.org/10.1007/s11356-012-0964-3> on accumulation of
% tetrazepam in _Mytilus galloprovincialis_. This set of files is a slight
% modification from the basic example file for the BYOM platform. Results
% (parameter estimates and CIs) are slightly different from those in the
% book chapter. I am not sure what has caused the discrepancy, but I think
% the accuracy of the ODE solver was not set high enough for the book's
% case study.
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

% observed internal concentrations (two treatments) over time (days)
% data extracted with PlotReader
DATA_extracted = [1 2.3	14.5
    0	0	0
    0	0	0
    0	0	0
    1	28.07029377	421.0526453
    1	49.66264269	598.110707
    1	79.89209129	883.1309136
    3	53.98111247	852.9015451
    3	101.4846003	982.456199
    3	120.9177944	NaN
    7	153.3063979	1209.176823
    7	194.332101	1543.859653
    7	228.8800194	1563.292847
    8	58.29974237	654.2510544
    8	71.25515173	766.5317491
    8	86.36987603	930.6343214
    10	15.1147243	582.9959827
    10	30.2294486	129.554734
    10	53.98111247	NaN
    14	0	174.8989069
    14	12.95556947	82.05140624
    14  NaN     NaN];

% data obtained from original author (only high treatment)
DATA_raw = [1 14.5
    0	0
    0	0
    0	0
    1	425
    1	598
    1	883
    3	853
    3	980
    3	856
    7	1208
    7	1548
    7	1562
    8	660
    8	933
    8	771
    10	586
    10	586
    10	135
    14	89
    14	82
    14	183];

DATA{1} = [DATA_extracted(:,[1 2]) DATA_raw(:,2)];

% if weight factors are not specified, ones are assumed in start_calc.m

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios.

X0mat = [2.3 14.5    % the scenarios (here exposure concentrations) 
         0    0  ];  % initial values state 2 (internal concentrations)

%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. 

% syntax: par.name = [startvalue fit(0/1) minval maxval];
par.ku    = [10   1 0 100];  % uptake rate constant, L/kg/d
par.ke    = [0.2  1 0 100];  % elimination rate constant, d-1

%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% used, based on the data set

glo.t   = linspace(0,16,100); % time vector for the model curves in days

% specify the y-axis labels for each state variable
glo.ylab{1} = ['internal concentration (',char(181),'g/kg dwt)'];
% specify the x-axis label (same for all states)
glo.xlab    = 'time (days)';
glo.leglab1 = 'Conc. '; % legend label before the 'scenario' number
glo.leglab2 = [char(181),'g/L']; % legend label after the 'scenario' number

prelim_checks % script to perform some preliminary checks and set things up
% Note: prelim_checks also fills all the options (opt_...) with defauls, so
% modify options after this call, if needed.

%% Calculations and plotting
% Here, the functions are called that will do the calculation and the
% plotting. 
% 
% Options for the plotting can be set using opt_plot (see prelim_checks.m).
% Options for the optimsation routine can be set using opt_optim. Options
% for the ODE solver are part of the global glo. 

glo.eventson   = 0; % events function on (1) or off (0)
glo.useode     = 1; % calculate model using ODE solver (1) or analytical solution (0)
opt_optim.fit  = 1; % fit the parameters (1), or don't (0)
opt_optim.it   = 1; % show iterations of the optimisation (1, default) or not (0)
opt_plot.bw    = 1; % if set to 1, plots in black and white with different plot symbols
opt_plot.annot = 0; % annotations in multiplot for fits: 1) box with parameter estimates 2) single legend

% optimise and plot (fitted parameters in par_out)
par_out = calc_optim(par,opt_optim); % start the optimisation
calc_and_plot(par_out,opt_plot); % calculate model lines and plot them

%% Profiling the likelihood
% By profiling you make robust confidence intervals for one or more of your
% parameters. Use the names of the parameters as they occurs in your
% parameter structure _par_ above. This can be a single string (e.g.,
% 'kd'), a cell array of strings (e.g., {'kd','ke'}), or 'all' to profile
% all fitted parameters. This example produces a profile for two parameters
% and provides the 95% confidence interval (on screen and indicated by the
% horizontal broken line in the plot).
% 
% *Note: for more post-calculations, see <byom_bioconc_extra.html
% byom_bioconc_extra.m>*
%
% Options for profiling can be set using opt_prof (see prelim_checks.m).

opt_prof.detail   = 2; % detailed (1) or a coarse (2) calculation
opt_prof.subopt   = 0; % number of sub-optimisations to perform to increase robustness

calc_proflik(par_out,{'all'},opt_prof,opt_optim);  % calculate a profile

%% Other files: derivatives
% To archive analyses, publishing them with Matlab is convenient. To keep
% track of what was done, the file derivatives.m can be included in the
% published result.
% 
% <include>derivatives.m</include>

%% Other files: call_deri
% To archive analyses, publishing them with Matlab is convenient. To keep
% track of what was done, the file call_deri.m can be included in the
% published result.
%
% <include>call_deri.m</include>