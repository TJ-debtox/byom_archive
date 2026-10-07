%% BYOM, byom_debtox_capitella_marecotox.m
% This is the file to reproduce the case study in Section 3.6.5 of the
% book (Fig. 3.8).
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
% *The model:* DEBtox following Jager & Zimmer (2012),
% <http://dx.doi.org/10.1016/j.ecolmodel.2011.11.012>. This set of files is
% taken from the BYOM 'DEBtox_classic' package that can be downloaded from
% <http://www.debtox.info/byom.html>.
%
% *This script:* This data set is from Hansen et al. (1999),
% <http://dx.doi.org/10.1890/1051-0761(1999)009%5B0482:EONNOL%5D2.0.CO;2>.
% The study is on _Capitella teleta_ exposed to nonylphenol in sediment.
% Analysed with full DEB earlier by Jager and Selck (2011),
% <http://dx.doi.org/10.1016/j.seares.2011.04.003>. Simultaneous fit on
% growth and reproduction.
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

% scaled damage
DATA{1} = [0];
        
% scaled reserve density
DATA{2} = [0];   

% body volume (mm3), first column in days
DATA{3} = [0.5 0	14	52	174
    0	0.0041	0.0041	0.0041	0.0041
    14	0.80	1.17	1.17	0.80
    21	3.81	4.73	3.44	2.21
    25	3.99	6.33	4.30	2.83
    32	8.48	12.16	7.74	4.36
    39	15.17	18.92	11.61	7.37
    46	15.42	19.53	14.19	10.26
    53	14.74	19.29	14.25	11.12
    60	16.83	19.16	16.89	13.02
    66	16.89	21.01	16.58	13.88
    72	18.24	22.17	16.65	16.09
    78	16.77	20.15	18.92	18.49];

DATA{3}(2:end,2:end) = DATA{3}(2:end,2:end).^(1/3); % convert volume to vol. length

% cumulative reproduction, first column in days
DATA{4} = [0.5 0	14	52	174
    14	0	0	0	0
    21	0	0	0	0
    25	0	0	0	0
    32	0	0	0	0
    39	493	604	343	0
    46	1007	1234	775	212
    53	1461	1768	1145	427
    60	1936	2284	1599	655
    66	2400	2914	2065	845
    72	2868	3441	2464	1033
    78	3256	3901	2876	1199];

% if weight factors are not specified, ones are assumed in start_calc.m

% Put the position of the various states in globals, to make sure correct
% one is selected for extra things (e.g., using plot_tktd for nicer plots);
% note that this must match the positions in the state vector in call_deri
% and derivatives.
glo.locD = 1; % location of scaled conc./damage in the state variable list
glo.locL = 3; % location of body size in the state variable list
glo.locR = 4; % location of cumulative reproduction in the state variable list
glo.locS = []; % location of survival probability in the state variable list

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios.

X0mat = [0 14 52 174   % the scenarios (here nominal concentrations) 
         0 0   0   0   % initial values state 1 
         1 1   1   1   % initial values state 2
         0 0   0   0   % initial values state 3 (overwritten by L0)
         0 0   0   0]; % initial values state 4
        
X0mat(:,[1]) = []; % remove the control treatment for 'hormesis'

%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. 
  
% global parameters as part of the structure glo
glo.moa = 6; % mode of action (see derivatives)

% syntax: par.name = [startvalue fit(0/1) minval maxval optional:log/normal scale (0/1)];
par.ke    = [0.4305 1   0    10 1]; % elimination rate constant (d-1)
par.c0    = [14     1  14   1e6 1]; % no-effect concentration (mg/kg)
% note that the lower boundary for c0 is set to the external concentration
% in the NEW control (first treatment)
par.cT    = [241.5  1  0    1e6 1]; % tolerance concentration (mg/kg)
par.c0s   = [1e6    0  0    1e6 1]; % no-effect concentration survival (uM)
par.b     = [0      0  0    1e6 1]; % killing rate (1/(d uM))

par.L0   = [0.16    0 1e-3 100 1]; % initial body length (mm)
par.Lp   = [2.051   1 1e-3 100 1]; % length at puberty (mm)
par.Lm   = [2.807   1 1e-3 100 1]; % maximum length (mm)
par.rB   = [0.0529  1    0 100 1]; % von Bertalanffy growth rate (d-1)
par.Rm   = [123.4   1    0 1e6 1]; % maximum reproduction rate (eggs/d)
par.g    = [10      0    0 100 1]; % energy-investment ratio (-)
par.f0   = [0.97    0    0   2 1]; % scaled functional response in control (-)
par.f    = [1       0    0   2 1]; % scaled functional response (-)
par.Tlag = [5.654   1    0  30 1]; % lag time for start development (d)

%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% used, based on the data set

% specify what to calculate and what to plot
glo.t = linspace(0,80,100); % time vector for the model curves in days

% specify the y-axis labels for each state variable
glo.ylab{1} = 'scaled damage (mg/kg)';
glo.ylab{2} = 'scaled reserve density (-)';
glo.ylab{3} = 'volumetric body length (mm)';
glo.ylab{4} = 'cumululative reproduction (eggs)';

% specify the x-axis label (same for all states)
glo.xlab    = 'time (days)';
glo.leglab1 = ''; % legend label before the 'scenario' number
glo.leglab2 = 'mg/kg'; % legend label after the 'scenario' number

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

glo.stiff = [0 3]; % ODE solver 0) ode45 (standard), 1) ode113 (moderately stiff), 2) ode15s (stiff)
% Second argument is for normally tight (1), tighter (2), or very tight (3)
% tolerances. Use 1 for quick analyses, but check with 3 to see if there is
% a difference! Especially for time-varying exposure, there can be large
% differences between the settings!
glo.break_time = 0; % break time vector up for ODE solver (1) or don't (0)
% Note: breaking the time vector is a good idea when the exposure scenario
% contains discontinuities. Don't use for continuous splines (type 1) as it
% will be much slower. For FOCUS scenarios (high time resolution), breaking
% up is not efficient and does not appear to be necessary.

opt_optim.fit    = 1; % fit the parameters (1), or don't (0)
opt_optim.it     = 1; % show iterations of the optimisation (1, default) or not (0)
opt_plot.bw      = 1; % if set to 1, plots in black and white with different plot symbols
opt_plot.annot   = 2; % annotations in multiplot for fits: 1) box with parameter estimates 2) single legend
opt_plot.statsup = [2]; % vector with states to suppress in plotting fits

% optimise and plot (fitted parameters in par_out)
par_out = calc_optim(par,opt_optim); % start the optimisation
calc_and_plot(par_out,opt_plot); % calculate model lines and plot them

% these plots are more readable, especially when plotting CIs
opt_tktd.repls = 0; % plot individual replicates (1) or means (0)
plot_tktd(par_out,opt_tktd,[]);
% Leave the options opt_conf empty to suppress all CIs for these plots.
% Plotting CIs requires a sample as saved by the various methods available
% in BYOM (see examples directory).

%% Profiling the likelihood
% By profiling you make robust confidence intervals for one or more of your
% parameters. Use the names of the parameters as they occurs in your
% parameter structure _par_ above. This can be a single string (e.g.,
% 'kd'), a cell array of strings (e.g., {'ke','c0'}), or 'all' to profile
% all fitted parameters. If you have used the parameter-space explorer for
% fitting, you already have (even more robust) profiles!
%
% Options for the profiling can be set using opt_prof (see prelim_checks).
% For this demo, the level of detail of the profiling is set to 'coarse'
% and no sub-optimisations are used. However, consider these options when
% working on your own data (e.g., set opt_prof.subopt=10).

opt_prof.detail   = 2; % detailed (1) or a coarse (2) calculation
opt_prof.subopt   = 0; % number of sub-optimisations to perform to increase robustness
opt_prof.brkprof  = 2; % when a better optimum is located, stop (1) or automatically refit (2)

% % UNCOMMENT FOLLOWING LINE(S) TO CALCULATE
% calc_proflik(par_out,{'ke','c0'},opt_prof,opt_optim);  % calculate a profile
% % Enter single parameter names, a cell array of names, or 'all' to profile
% % all fitted parameters.

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