%% BYOM, byom_snail_fit.m, DEBkiss with Lymnaea 
%
% * Author: Tjalling Jager 
% * Date: November 2021
% * Web support: <http://www.debtox.info/byom.html>
% * Back to index <walkthrough_debkiss_paper.html>
%
% BYOM is a General framework for simulating model systems in terms of
% ordinary differential equations (ODEs). The model itself needs to be
% specified in <derivatives.html derivatives.m>, and <call_deri.html
% call_deri.m> may need to be modified to the particular problem as well.
% The files in the engine directory are needed for fitting and plotting.
% Results are shown on screen but also saved to a log file (results.out).
%
% *The model:* DEBkiss following Jager et al (2013),
% <http://dx.doi.org/10.1016/j.jtbi.2013.03.011>.
%
% *This script:* The case study in the original DEBkiss paper for the pond
% snail _Lymnaea stagnalis_ at three food levels. The data set was also
% used in Zimmer et al (2012),
% <http://dx.doi.org/10.1007/s10646-012-0973-5>.
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

% Shell length in mm over time in three feeding regimes
DATA{1} = [0.5        1            2             3
            0       12.918       12.799       12.618
           14        18.65       16.717       14.479
           28       21.549       19.385       17.019
           42       23.759       21.357       19.247
           56       25.669       23.372       21.095
           70       27.146       24.552       22.234
           84       27.953       25.514        22.95
           98       28.189       25.928       23.238
          112        28.27       26.317       23.709
          128       29.242       26.614       24.057
          140        29.53       26.819       24.153];
        
% weight factors (number of replicates per observation)
W{1} = [30    30    30
    30    30    29
    29    30    29
    28    29    28
    28    29    28
    26    29    28
    26    28    28
    26    28    28
    25    28    26
    23    28    26
    21    28    26];

% Cumulative reproduction over time in three feeding regimes; reduced data
% to match the frequency of the size measurements (see DEBkiss paper)
DATA{2}= [0.5            1            2            3
           35            0            0            0
           49           76         17.3          0.8
           63        219.8         79.3          7.9
           77        389.7        153.4         42.6
           91        534.6        238.6           82
          105        676.8        308.6        135.4
          119        858.9        413.9        187.7
          133       1023.1        514.1        249.6];

% weight factors (number of replicates per observation)
W{2} = [28    29    28
    28    29    27
    28    29    27
    26    29    27
    26    28    27
    26    28    27
    24    28    26
    23    28    26];

DATA{3}=0; % no data for the egg buffer (not used in this script)

% if weight factors are not specified, ones are assumed in start_calc.m

% Notes: In this example, the time vectors are not the same for both data 
% sets. This is no problem for the analysis. 

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios.

X0mat = [   1    2    3    % the scenarios (here feeding level) 
         12.8 12.8 12.8    % initial shell length in mm
            0    0    0    % initial cumul repro
            0    0    0 ]; % initial weight of buffer in egg

% Put the position of the various states in globals, to make sure correct
% one is selected for extra things (e.g., to prevent shrinking in
% call_deri, and for using plot_tktd).
glo.locL = 1; % location of body size in the state variable list
glo.locR = 2; % location of cumulative reproduction in the state variable list

%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. 
% Global parameters as part of the structure glo

glo.delM  = 0.401; % shape corrector (used in call_deri.m)
glo.dV    = 0.1;   % dry weight density (used in call_deri.m)
glo.len   = 2;     % switch to fit physical length (0=off, 1=on, 2=on and no shrinking) (used in call_deri.m)
glo.Tlag  = 0;     % delay for start development (only used for embryo)
glo.mat   = 0;     % include maturity maint. (0=off, 1=include)
% Note: if glo.len > 0 than the initial state for size in X0mat is length
% too (see call_deri.m)

% syntax: par.name = [startvalue fit(0/1) minval maxval optional:log/normal scale (0/1)];
par.sJAm = [0.111   1 0 1e6]; % specific assimilation rate 
par.sJM  = [0.00795 1 0 1e6]; % specific maintenance costs 
par.WB0  = [0.15    0 0 1e6]; % initial dry weight of egg
par.WVp  = [70.3    1 0 1e6]; % body mass at puberty
par.yAV  = [0.8     0 0 1];   % yield of assimilates on volume (starvation)
par.yBA  = [0.95    0 0 1];   % yield of egg buffer on assimilates
par.yVA  = [0.8     0 0 1];   % yield of structure on assimilates (growth)
par.kap  = [0.893   1 0 1];   % allocation fraction to soma
par.fB   = [1       0 0 2];   % scaled food level, embryo (not used here)
par.f1   = [1       0 0 2];   % scaled food level, ad libitum
par.f2   = [0.891   1 0 2];   % scaled food level, regime 2
par.f3   = [0.797   1 0 2];   % scaled food level, regime 3
par.WVf  = [0       0 0 1e6]; % half-saturation body weight

% % Adding the male function
% par.sJAm = [0.119 1 0 1e6]; % specific assimilation rate 
% par.yBA  = [0.55  0 0 1];   % yield of egg buffer on assimilates (lowered!)
% par.kap  = [0.828 1 0 1];   % allocation fraction to soma

% % Adding maturity maintenance yields a different fit
% glo.mat  = 1;     % include maturity maint. (0=off, 1=include)
% par.sJAm = [0.12   1 0 1e6]; % specific assimilation rate 
% par.sJM  = [0.0078 1 0 1e6]; % specific maintenance costs 
% par.WVp  = [67  1 0 1e6];    % body mass at puberty
% par.kap  = [0.776  1 0 1];   % allocation fraction to soma
% par.f2   = [0.904  1 0 2];   % scaled food level, regime 2
% par.f3   = [0.831  1 0 2];   % scaled food level, regime 3
% par.yBA  = [0.95   0 0 1];   % yield of egg buffer on assimilates

% % Adding initial food limitation
% glo.mat   = 1;     % include maturity maint. (0=off, 1=include)
% par.sJAm = [0.18  1 0 1e6]; % specific assimilation rate 
% par.sJM  = [0.011 1 0 1e6]; % specific maintenance costs 
% par.WVp  = [60.3  1 0 1e6]; % body mass at puberty
% par.kap  = [0.77  1 0 1];   % allocation fraction to soma
% par.f2   = [0.904 1 0 2];   % scaled food level, regime 2
% par.f3   = [0.831 1 0 2];   % scaled food level, regime 3
% par.yBA  = [0.55  0 0 1];   % yield of egg buffer on assimilates
% par.WVf  = [5     1 0 1e6]; % half-saturation body weight

%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% used, based on the data set

glo.t   = linspace(0,145,100); % time vector for the model simulation in days

% specify the y-axis labels for each state variable
glo.ylab{1} = 'shell length (mm)';
glo.ylab{2} = 'cumulative reproduction (eggs)';
glo.ylab{3} = 'egg buffer (mg)';

% specify the x-axis label (same for all states)
glo.xlab    = 'time (days)';
glo.leglab1 = 'feeding lvl. '; % legend label before the 'scenario' number
glo.leglab2 = ''; % legend label after the 'scenario' number

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
% 
% For the demo, the iterations are turned off (opt_optim.it = 0).

glo.stiff = [0 3]; % ODE solver 0) ode45 (standard), 1) ode113 (moderately stiff), 2) ode15s (stiff)
% Second argument is for normally tight (1), tighter (2), or very tight (3)
% tolerances. Use 1 for quick analyses, but check with 3 to see if there is
% a difference!

opt_optim.fit    = 1; % fit the parameters (1), or don't (0)
opt_optim.it     = 1; % show iterations of the optimisation (1, default) or not (0)
opt_plot.bw      = 0; % if set to 1, plots in black and white with different plot symbols
opt_plot.annot   = 2; % annotations in multiplot for fits: 1) box with parameter estimates 2) single legend
opt_plot.statsup = [3]; % vector with states to suppress in plotting fits

% optimise and plot (fitted parameters in par_out)
par_out = calc_optim(par,opt_optim); % start the optimisation
calc_and_plot(par_out,opt_plot); % calculate model lines and plot them

%% Profiling the likelihood
% By profiling you make robust confidence intervals for one or more of your
% parameters. Use the names of the parameters as they occurs in your
% parameter structure _par_ above. This can be a single string (e.g.,
% 'kd'), a cell array of strings (e.g., {'rB','Lm'}), or 'all' to profile
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
% calc_proflik(par_out,{'sJAm','sJM'},opt_prof,opt_optim);  % calculate a profile
% % Enter single parameter names, a cell array of names, or 'all' to profile
% % all fitted parameters.

%% Population growth rate
% The function calc_pop allows for a calculation of the intrinsic rate of
% population increase (Euler-Lotka equation). It calculates the growth rate
% at three food levels (set in the engine file calc_pop.m). Note that two
% plots are made: the absolute growth rate and the relative one.

% Population calculations can be performed with CIs
opt_conf.type    = 0; % make intervals from 1) slice sampler, 2)likelihood region, 3) parspace explorer
opt_conf.lim_set = 0; % for lik-region sample: use limited set of n_lim points (1) or outer hull (2) to create CIs

Tpop = linspace(0,300,100); % time vector for the population calculations
% Cpop = linspace(0,1,50);   % concentration range for population calculations
Cpop = -1; % use only scenarios as given in X0mat
Th   = 10; % time from fresh egg to t=0 in time vector (e.g., hatching time)
glo.locS = []; % no survival data (calc_pop needs to know that)

% UNCOMMENT FOLLOWING LINE(S) TO CALCULATE
calc_pop(par_out,X0mat,Tpop,Cpop,Th,opt_conf,opt_pop)

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
  