%% BYOM, byom_gutsred_timevar_diazinon.m, reduced GUTS model
%
% * Author: Tjalling Jager 
% * Date: November 2021
% * Web support: <http://www.debtox.info/byom.html>
% * Back to index <walkthrough_guts.html>
%
% BYOM is a General framework for simulating model systems. The files in
% this directory use the ODE solver only, and therefore <simplefun.html
% simplefun.m> will be missing.
%
% *The model:* fitting survival data with the
% <http://www.debtox.info/about_guts.html GUTS> special cases based on the
% reduced model (TK and damage dynamics lumped): SD, IT and mixed (or GUTS
% proper). 
% 
% *This script:* dataset from Ashauer et al. (2010),
% <http://dx.doi.org/10.1021/es903478b>. _Gammarus pulex_ exposed to pulses
% of diazinon. One scenario from this data set (two pulses) was used in
% Jager (2014), <http://dx.doi.org/10.1007/s10646-013-1149-7>. This script
% demonstrates two ways to define time-varying exposure scenarios: linear
% interpolation or a series of block pulses. The files in this directory
% use the ODE solver only!
% 
%  Copyright (c) 2012-2023, Tjalling Jager.
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

pathdefine(1) % set path to the BYOM/engine directory (option 1 uses parallel toolbox)
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

% observed number of survivors, time in days, conc in nM
DATA{1} = [-1     0     1     2     3
     0    60    70    70    70
     1    55    66    65    65
     2    53    61    59    59
     3    51    55    56    55
     4    51    31    54    53
     5    51    31    50    51
     6    51    29    47    48
     7    48    26    46    46
     8    48    24    46    46
     9    46    22    40    46
    10    45    21    23    44
    11    44    19    22    41
    12    42    17    22    40
    13    41    14    21    40
    14    41    14    18    40
    15    40    13    17    39
    16    38    11    17    38
    17    38    11    13    36
    18    37    10    13    33
    19    37     9    13    28
    20    37     8    11    24
    21    37     8    11    23
    22    36     8    11    19];

% In this data set, exposure was time-varying and reported as a series of
% concentrations over time. Here, the scenario is used as a linear forcing
% series (which has an analytical solution, and is thus much faster than
% the ODE version). Double time entries are used, which is more efficient,
% and probably more accurate.
Cw1 = [0    1
    0   102.65
    1	97.59
    1	0
    3	0
    3	103.88
    4	98.19
    4	0
    22	0];

Cw2 = [0    2
    0   100.78
    1	106.32
    1	0
    8   0
    8	103.56
    9   95.82
    9	0
    22	0];

Cw3 = [0    3
    0   100.6
    1	94.61
    1	0
    16  0
    16	100.58
    17  96.51
    17	0
    22	9.85];

make_scen(4,Cw1,Cw2,Cw3); % prepare as linear-forcing function interpolation (set glo.use_ode = 0)

% Note: it is also possible to use smooth splines (type 1), but you might
% as well set glo.break_time = 0, since it won't do much.

% Create a table with nicer custom labels for the legends
Scenario = [0;1;2;3];
Label = {'control';'2-d interval'; '7-d interval';'15-d interval'};
glo.LabelTable = table(Scenario,Label); % create a Matlab table for the labels

% internal concentrations
DATA{2} = 0; % no data for scaled body residues

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios.

X0mat(1,:) = [0 1 2 3]; % scenarios (concentrations or identifiers)
X0mat(2,:) = 1;         % initial survival probability
X0mat(3,:) = 0;         % initial internal concentration

%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. 

% global parameters for GUTS purposes
glo.sel  = 1; % select death mechanism: 1) SD 2) IT 3) mixed
glo.locD = 2; % location of scaled damage in the state variable list
glo.locS = 1; % location of survival probability in the state variable list
 
% Best values SD
% syntax: par.name = [startvalue fit(0/1) minval maxval optional:log/normal scale (0/1)];
par.kd    = [0.083 1 1e-3  10 0]; % dominant rate constant (d-1)
par.mw    = [4.8   1    0 1e6 1]; % median threshold for survival (nM)
par.hb    = [0.027 1    0 1e6 1]; % background hazard rate (1/d)
par.bw    = [0.023 1 1e-6 1e3 0]; % killing rate (nM-1 d-1) (SD and mixed)
par.Fs    = [2     1    1 100 1]; % fraction spread of NEC distribution (-) (IT and mixed)

% % Best values IT
% % syntax: par.name = [startvalue fit(0/1) minval maxval optional:log/normal scale (0/1)];
% par.kd    = [0.01 1 1e-3  10 0]; % dominant rate constant (d-1)
% par.mw    = [2.3   1    0 1e6 1]; % median threshold for survival (nM)
% par.hb    = [0.051 1    0 1e6 1]; % background hazard rate (1/d)
% par.bw    = [0.023 1    0 1e3 0]; % killing rate (nM-1 d-1) (SD and mixed)
% par.Fs    = [1.1     1    1 100 1]; % fraction spread of NEC distribution (-) (IT and mixed)

switch glo.sel % make sure that right parameters are fitted
    case 1 % SD
        par.Fs(2) = 0; % never fit the NEC spread!
    case 2 % IT
        par.bw(2) = 0; % never fit the killing rate!
    case 3 % mixed
        % do nothing: fit all parameters
end

%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% constructed, based on the data set.

% specify the y-axis labels for each state variable
glo.ylab{1} = 'survival probability';
glo.ylab{2} = 'scaled damage (nM)';
% specify the x-axis label (same for all states)
glo.xlab    = 'time (days)';
glo.leglab1 = 'scenario '; % legend label before the 'scenario' number
glo.leglab2 = 'nM'; % legend label after the 'scenario' number

prelim_checks % script to perform some preliminary checks and set things up
% Note: prelim_checks also fills all the options (opt_...) with defauls, so
% modify options after this call, if needed.

%% Calculations and plotting
% Here, the function is called that will do the calculation and the plotting.
% Options for the plotting can be set using opt_plot (see prelim_checks.m).
% Options for the optimsation routine can be set using opt_optim. Options
% for the ODE solver are part of the global glo. 
% 
% NOTE: for this package, the options useode and eventson in glo will not
% be functional: the ODE solver is always used, and the events function as
% well.

glo.stiff    = [0 1]; % ODE solver 0) ode45 (standard), 1) ode113 (moderately stiff), 2) ode15s (stiff)
% Seems that standard solver performs good enough, with tight tolerances.
% Second argument is for normally tight (1), tighter (2), or very tight (3)
% tolerances.
glo.break_time = 1; % break time vector up for ODE solver (1) or don't (0)

opt_plot.annot = 2; % annotations in multiplot for fits: 1) box with parameter estimates 2) single legend
opt_optim.fit  = 1; % fit the parameters (1), or don't (0)
opt_optim.it   = 1; % show iterations of the optimisation (1, default) or not (0)

% optimise and plot (fitted parameters in par_out)
par_out = calc_optim(par,opt_optim); % start the optimisation
calc_and_plot(par_out,opt_plot); % calculate model lines and plot them

opt_tktd.obspred = 0; % plot predicted-observed plots (1) or not (0)
plot_tktd(par_out,opt_tktd,[]); % leaving opt_conf empty suppresses all CIs for these plots

% Below, several options for post analyses are provided. These will be very
% slow when using the ODE solver! The parallel toolbox helps here.

return % stop here, you can run the next sections afterwards if you like

%% Profiling the likelihood
% By profiling you make robust confidence intervals for one or more of your
% parameters. 

opt_prof.detail   = 2; % detailed (1) or a coarse (2) calculation
opt_prof.subopt   = 10; % number of sub-optimisations to perform to increase robustness
opt_prof.brkprof  = 2; % when a better optimum is located, stop (1) or automatically refit (2)

% % UNCOMMENT FOLLOWING LINE(S) TO CALCULATE
% calc_proflik(par_out,'all',opt_prof,opt_optim);  % calculate a profile
% % Enter single parameter names, a cell array of names, or 'all' to profile
% % all fitted parameters.

%% Likelihood region
% Make intervals on model predictions by using a sample of parameter sets
% taken from the joint likelihood-based conf. region. 

opt_prof.detail  = 2;  % detailed (1) or a coarse (2) calculation
opt_prof.subopt  = 20; % number of sub-optimisations to perform to increase robustness
opt_prof.brkprof = 2; % when a better optimum is located, stop (1) or automatically refit (2)
opt_likreg.burst = 5000; % number of random sets from parameter space evaluated every iteration
opt_likreg.skipprof = 0; % skip profiling step; use boundaries from saved likreg set (1) or profiling (2)

par_better = calc_likregion(par_out,500,opt_likreg,opt_prof,opt_optim); % second argument is target for number of samples in inner region (-1 to re-use saved sample from previous runs)
% Second entry is the number of accepted parameter sets to aim for. Use -1
% here to use a saved set.

if ~isempty(par_better) % if the profiling found a better optimum ...
    calc_and_plot(par_better,opt_plot); % calculate model lines and plot them
    par_out = par_better; % use the new parameter structure for further analyses below
end

%% Plot results with confidence intervals
% The following code can be used to make plots with confidence intervals.
% Options for confidence bounds on model curves can be set using opt_conf
% (see prelim_checks). The plot_tktd function makes multiplots for the
% survival data, which are more readable when plotting with various
% intervals.

opt_conf.type    = 2; % make intervals from 1) slice sampler, 2)likelihood region, 3) parspace explorer
opt_conf.lim_set = 2; % for lik-region sample: use limited set of n_lim points (1) or outer hull (2) to create CIs
opt_tktd.obspred = 1; % plot predicted-observed plots (1) or not (0)

plot_tktd(par_out,opt_tktd,opt_conf); 

%% Calculate LCx versus time
% Here, the LCx (by default the LC50) is calculated at several time points.
% LCx values are also printed on screen. If a sample from parameter space
% is available (e.g., from the slice sampler or the likelihood region), it
% can be used to calculate confidence bounds. 
% 
% Options for LCx (with confidence bounds) can be set using opt_ecx (see
% prelim_checks). Note that opt_conf.type=-1 skips CIs.

opt_conf.type    = -1; % make intervals from 1) slice sampler, 2)likelihood region, 3) parspace explorer
opt_conf.lim_set = 2; % for lik-region sample: use limited set of n_lim points (1) or outer hull (2) to create CIs

% This is the slower general method, which allows for changes in the model
% in derivatives.m. The fast method used in the 'reduced' folder makes use
% of the analytical solutions, so they are best restricted to the 'reduced'
% folder,
opt_ecx.Feff      = [0.50]; % effect levels (>0 en <1), x/100 in ECx
opt_ecx.notitle   = 1; % set to 1 to suppress titles above ECx plots
Tend = [2:0.2:3 3.5 4:8 10 15 22]; % times at which to calculate LCx, relative to control

calc_ecx(par_out,Tend,opt_ecx,opt_conf); % general method for ECx values

%% Calculate LPx for a monitoring profile
% We can use the function calc_lpx_lim to calculate an exposure
% multiplication factor: with which factor do we need to multiply an
% exposure scenario to obtain x% effect at the end of the scenario? Here,
% done for 10% effect at the end of the exposure scenario. Note that
% opt_conf.type = -1 skips CIs. Also not that profile_monit.txt needs to be
% in the current directory (or the search path), and needs to be a simple
% two-column scenario (time points and concentrations). It is also possible
% to use a scenario identifier here (see the diazinon case study in this
% directory).

opt_conf.type     = -1; % make intervals from 1) slice sampler, 2)likelihood region shooting, 3) parspace explorer
opt_conf.lim_set  = 2; % for lik-region sample: use limited set of n_lim points (1) or outer hull (2) to create CIs
opt_ecx.Feff      = [0.10]; % effect levels (>0 en <1), x/100 in ECx
opt_ecx.notitle   = 1; % set to 1 to suppress titles above ECx plots

glo.break_time = 0; % break time vector up for ODE solver (1) or don't (0)
% For relatively short profiles, there is still some gain in calculation
% time from breaking the time vector. For 485-day FOCUS profile, there is
% hardly any gain. However, breaking up should have higher precision.

% This is the slower general method, which allows for changes in the model
% in derivatives.m. The fast method used in the 'reduced' folder makes use
% of the analytical solutions, so they are best restricted to the 'reduced'
% folder.

calc_epx(par_out,'profile_monit.txt',[],opt_ecx,opt_conf,opt_tktd); % general method for EPx values
% calc_epx(par_out,'profile_focus.txt',[],opt_ecx,opt_conf,opt_tktd); % general method for EPx values
