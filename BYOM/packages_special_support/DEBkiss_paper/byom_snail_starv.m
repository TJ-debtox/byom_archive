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
% snail _Lymnaea stagnalis_ under starvation. There is no fitting in this
% script. Data from the MD snails of Zonneveld and Kooijman (1989).
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

DATA{2}=0; % no data for reproduction
DATA{3}=0; % no data for egg buffer

% body weight on starvation (mg dwt)
DATA{1} = [1 2
    0	185.5537113
    0	168.3152507
    0	154.6254851
    0	171.3550357
    7	163.6021266
    7	171.9668996
    7	178.0503701
    7	184.6407979
    7	149.6588789
    14	167.5072485
    14	175.3650664
    14	164.717036
    14	145.450755
    14	150.0133601
    22	158.9684937
    22	162.2637065
    22	128.2956999
    22	107.001594
    22	126.7728798
    28	118.5345235
    28	146.1636181
    28	172.2718519
    28	138.0503724
    35	140.4346243
    35	129.0261658
    35	79.08903684
    35	98.09989452
    35	66.66666812
    44	88.5276719
    44	97.39940913
    44	111.3406975
    44	105.00179
    44	108.0415749
    51	91.41497836
    51	98.25888603
    51	61.24915341
    51	92.68041588
    51	68.85152513
    58	68.19600142
    58	73.51903465
    58	103.682923
    58	59.32232746
    58	65.6592618
    65	53.0902793
    65	74.12699336
    65	76.66177803
    65	77.92916126
    65	81.47785615
    65	69.05547832];

W0 = mean(DATA{1}(2:5,2)); % initial average dwt

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios.

X0mat = [   1  2 % the scenarios (here two embryo feeding levels) 
            W0 W0   % initial structural dry weight (mg) (note: glo.len = 0)
            0  0    % initial cumul repro
            0  0 ]; % initial weight of buffer in egg (mg)

% Put the position of the various states in globals, to make sure correct
% one is selected for extra things (e.g., to prevent shrinking in
% call_deri, and for using plot_tktd).
glo.locL = 1; % location of body size in the state variable list
glo.locR = 2; % location of cumulative reproduction in the state variable list

%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. 
% global parameters as part of the structure glo

glo.delM  = 0.401; % shape corrector (used in call_deri.m)
glo.dV    = 0.1;   % dry weight density (used in call_deri.m)
glo.len   = 0;     % switch to fit physical length (0=off, 1=on, 2=on and no shrinking) (used in call_deri.m)
glo.Tlag  = 0;     % Delay for start development (only used for embryo)
glo.mat   = 0;     % include maturity maint. (0=off, 1=include)
% Note: if glo.len > 0 than the initial state for size in X0mat is length
% too (see call_deri.m)

% syntax: par.name = [startvalue fit(0/1) minval maxval optional:log/normal scale (0/1)];
par.sJAm = [0.119   0 0 1e6]; % specific assimilation rate 
par.sJM  = [0.00795 0 0 1e6]; % specific maintenance costs 
par.WB0  = [0.15    0 0 1e6]; % initial dry weight of egg
par.WVp  = [70.3    0 0 1e6]; % body mass at puberty
par.yAV  = [0.8     0 0 1];   % yield of assimilates on volume (starvation)
par.yBA  = [0.55    0 0 1];   % yield of egg buffer on assimilates
par.yVA  = [0.8     0 0 1];   % yield of structure on assimilates (growth)
par.kap  = [0.828   0 0 1];   % allocation fraction to soma
par.fB   = [0.5     0 0 2];   % scaled food level, embryo
par.f1   = [0       0 0 2];   % scaled food level, ad libitum
par.f2   = [0       0 0 2];   % scaled food level, regime 2
par.f3   = [0       0 0 2];   % scaled food level, regime 3
par.WVf  = [0       0 0 1e6]; % half-saturation body weight

%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% used, based on the data set

glo.t = linspace(0,70,100); % time vector for the model simulation in days

% specify the y-axis labels for each state variable
glo.ylab{1} = 'body dry weight (mg)';
glo.ylab{2} = 'cumulative reproduction (eggs)';
glo.ylab{3} = 'egg buffer (mg)';

% specify the x-axis label (same for all states)
glo.xlab    = 'time (days)';
glo.leglab1 = 'scen. '; % legend label before the 'scenario' number
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

glo.stiff = [0 3]; % ODE solver 0) ode45 (standard), 1) ode113 (moderately stiff), 2) ode15s (stiff)
% Second argument is for normally tight (1), tighter (2), or very tight (3)
% tolerances. Use 1 for quick analyses, but check with 3 to see if there is
% a difference!

opt_optim.fit    = 0; % fit the parameters (1), or don't (0)
opt_optim.it     = 0; % show iterations of the optimisation (1, default) or not (0)
opt_plot.bw      = 0; % if set to 1, plots in black and white with different plot symbols
opt_plot.annot   = 0; % annotations in multiplot for fits: 1) box with parameter estimates 2) single legend
opt_plot.statsup = [2 3]; % vector with states to suppress in plotting fits

% optimise and plot (fitted parameters in par_out)
par_out = calc_optim(par,opt_optim); % start the optimisation
calc_and_plot(par_out,opt_plot); % calculate model lines and plot them

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

