%% BYOM, simple-compound model: byom_debtox_daphnia.m
%
% *Table of contents*

%% About
% * Author     : Tjalling Jager
% * Date       : September 2022
% * Web support: <http://www.debtox.info/byom.html>
% * Back to index <walkthrough_debtox2019.html>
%
% BYOM is a General framework for simulating model systems in terms of
% ordinary differential equations (ODEs). The model itself needs to be
% specified in <derivatives.html derivatives.m>, and <call_deri.html
% call_deri.m> may need to be modified to the particular problem as well.
% The files in the engine directory are needed for fitting and plotting.
% Results are shown on screen but also saved to a log file (results.out).
%
% *The model:* simple DEBtox model for toxicants, based on DEBkiss and
% formulated in compound parameters. The model includes flexible modules
% for toxicokinetics/damage dynamics and toxic effects. The DEBkiss e-book
% (see <http://www.debtox.info/book_debkiss.html>) provides a description
% of the model, as well as the publication of Jager in Ecological
% Modelling: <https://doi.org/10.1016/j.ecolmodel.2019.108904>.
%
% *This script:* The water flea _Daphnia magna_ exposed to fluoranthene;
% data set from Jager et al (2010),
% <http://dx.doi.org/10.1007/s10646-009-0417-z>. Published as case study in
% Jager & Zimmer (2012). Simultaneous fit on growth, reproduction and
% survival. Parameter estimates differ slightly from the ones in Jager &
% Zimmer due to several model changes and choices in the analysis.
% 
%  Copyright (c) 2012-2022, Tjalling Jager.
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

% scaled damage
DATA{1} = [0]; % there are never data for this state

% body length (mm) on each observation time (d), concentrations in uM
DATA{2} = [1 0 0 0.213	0.426	0.853
    0	0.88	0.88	0.88	0.88	0.88
    2	1.38	1.34	1.4     1.44	1.32
    4	1.96	1.72	1.9     1.88	1.74
    6	2.28	2.38	2.14	2.14	2.08
    8	2.34	2.52	2.56	2.36	2.16
    10	2.64	2.56	2.46	2.46	2.48
    12	2.66	2.56	2.6     2.58	2.7
    14	2.68	2.7     2.76	2.78	2.8
    16	2.88	2.78	2.82	NaN     NaN
    18	2.94	2.84	2.92	3.08	2.9
    20	3.06	3.02	3.02	3.22	2.76
    21	3.26	3.04	3.06	3       2.84];

% W{2} = 5 * ones(size(DATA{2})-1); % weights: 5 animals in each treatment

W{2} = 0;

% cumulative reproduction (nr. offspring per mother) (mm) on each
% observation time (d), concentrations in uM
DATA{3} = [1 0 0 	0.213	0.426	0.853
    0	0	0	0	0	0
    2	0	0	0	0	0
    4	0	0	0	0	0
    6	0	0	0	0	0
    8	1.2     2.1     0       0       0
    10	6.5     7.7     10.0	0.8     1.7
    12	17.9	25.4	16.3	2.0     1.7
    14	31.8	29.5	33.3	6.4     1.7
    16	50.6	48.8	56.2	9.4     1.7
    18	61.3	62.5	62.5	13.8	1.7
    20	81.0	76.4	64.1	16.4	1.7
    21	81.0	76.4	64.1	16.4	1.7];

% weights: number of individuals alive
W{3} = [10	10	10	10	10
	10	10	10	10	10
	10	10	10	10	10
	10	10	10	10	10
	10	10	10	10	10
	10	10	10	10	9
	10	10	10	10	7
	10	10	10	10	5
	10	10	10	10	5
	10	9	10	10	3
	10	9	10	9	3
	10	8	10	9	2]; 

% Global setting glo.Tbp is used to shift model output by three days. This
% means that we don't have to change the data set to account for
% brood-pouch delay.
glo.Tbp = 3; % shift model output for brood-pouch delay

% survivors on each observation time (d), concentrations in uM
DATA{4} = [-1	0	0	0.213	0.426	0.853
    0	10	10	10	10	10
    2	10	10	10	10	10
    4	10	10	10	10	10
    6	10	10	10	10	10
    8	10	10	10	10	10
    10	10	10	10	10	9
    12	10	10	10	10	7
    14	10	10	10	10	5
    16	10	10	10	10	5
    18	10	9	10	10	3
    20	10	9	10	9	3
    21	10	8	10	9	2];

% if weight factors are not specified, ones are assumed in start_calc.m

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios.

X0mat = [0 0.213 0.426 0.853]; % the scenarios (here nominal concentrations) 
X0mat(2,:) = 0; % initial values state 1 (scaled damage)
X0mat(3,:) = 0; % initial values state 2 (body length, initial value overwritten by L0)
X0mat(4,:) = 0; % initial values state 3 (cumulative reproduction)
X0mat(5,:) = 1; % initial values state 4 (survival probability)

% Put the position of the various states in globals, to make sure that the
% correct one is selected for extra things (e.g., for plotting in
% plot_tktd, in call_deri for accommodating 'no shrinking', for population
% growth rate).
glo.locD = 1; % location of scaled damage in the state variable list
glo.locL = 2; % location of body size in the state variable list
glo.locR = 3; % location of cumulative reproduction in the state variable list
glo.locS = 4; % location of survival probability in the state variable list

%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. 
  
% global parameters as part of the structure glo
glo.FBV = 0.02;    % dry weight of egg as fraction of structural body weight (-) (for losses with reproduction)
glo.KRV = 1;       % part. coeff. repro buffer and structure (kg/kg) (for losses with reproduction)
glo.kap = 0.8;     % approximation for kappa (for starvation response)
glo.yP  = 0.8*0.8; % product of yVA and yAV (for starvation response)
glo.len = 2;       % set switch to 2 to prevent shrinking in length (used in call_deri.m)

% select mode of action of toxicant as set of switches:
% [assimilation/feeding, maintenance costs (somatic and maturity), growth costs, repro costs] 
% glo.moa = [1 0 0 0 0]; % assimilation/feeding
% glo.moa = [0 1 0 0 0]; % costs for maintenance 
% glo.moa = [0 0 1 1 0]; % costs for growth and reproduction
glo.moa = [0 0 0 1 0]; % costs for reproduction
% glo.moa = [0 0 0 0 1]; % hazards for reproduction

% select which feedbacks to use on damage dynamics as set of switches:
% [surface:volume on uptake, surface:volume on elimination, growth dilution, losses with reproduction] 
% glo.feedb = [1 1 1 1]; % all feedbacks 
glo.feedb = [1 1 1 0]; % classic DEBtox (no losses with repro)
% glo.feedb = [0 0 1 0]; % damage that is diluted by growth
% glo.feedb = [0 0 0 0]; % damage that is not diluted by growth

% syntax: par.name = [startvalue fit(0/1) minval maxval optional:log/normal scale (0/1)];
par.L0   = [1     1 0 1e6 1]; % body length at start experiment (mm)
par.Lp   = [1.8   1 0 1e6 1]; % body length at puberty (mm)
par.Lm   = [3     1 0 1e6 1]; % maximum body length (mm)
par.rB   = [0.15  1 0 1e6 1]; % von Bertalanffy growth rate constant (1/d)
par.Rm   = [10    1 0 1e6 1]; % maximum reproduction rate (#/d)
par.f    = [1     0 0   2 1]; % scaled functional response
par.hb   = [0.005 1 1e-6 1e6 0]; % background hazard rate (d-1)
par.a    = [1     0 0.1  10  0]; % coefficient for Weibull backgound hazard (-)
% Note: it does not matter whether length measures are entered as actual
% length or as volumetric length (as long as the same measure is used
% consistently, and the organism is isomorphic).
% Note: hb is fitted on log-scale. This is especially helpful when hb is
% fitted along with the other parameters. The simplex fitting sometimes has
% trouble when one parameter is much smaller than others.

% extra parameters for special situations
par.Lf   = [0 0 0 1e6 1]; % body length at half-saturation feeding (mm)
par.Lj   = [0 0 0 1e6 1];  % body length at end acceleration (mm)
par.Tlag = [0 0 0  30 1]; % lag time for start development

ind_tox = length(fieldnames(par))+1; % index where tox parameters start

% the parameters below this line are all treated as toxicity parameters!
par.kd   = [0.076  1 0.01  10 0]; % dominant rate constant (d-1)
par.zb   = [0.096  1 0    1e6 1]; % effect threshold energy budget ([C])
par.bb   = [74     1 1e-6 1e6 0]; % effect strength energy-budget effect (1/[C])
par.zs   = [0.23   1 0    1e6 1]; % effect threshold survival ([C])
par.bs   = [0.66   1 1e-6 1e6 0]; % effect strength survival (1/([C] d))

%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% constructed, based on the data set.

% specify the y-axis labels for each state variable
glo.ylab{1} = ['scaled damage (',char(181),'M)'];
glo.ylab{2} = 'body length (mm)';
if isfield(glo,'Tbp') && glo.Tbp > 0
    glo.ylab{3} = ['cumul. repro. (shift ',num2str(glo.Tbp),'d)'];
else
    glo.ylab{3} = 'cumul. repro. (no shift)';
end
glo.ylab{4} = 'survival fraction (-)';

% specify the x-axis label (same for all states)
glo.xlab    = 'time (days)';
glo.leglab1 = 'conc. '; % legend label before the 'scenario' number
glo.leglab2 = [char(181),'M']; % legend label after the 'scenario' number

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
opt_plot.bw      = 0; % if set to 1, plots in black and white with different plot symbols
opt_plot.annot   = 2; % annotations in multiplot for fits: 1) box with parameter estimates 2) single legend
opt_plot.repls   = 0; % set to 1 to plot replicates, 0 to plot mean responses

% Note: using opt_optim.type = 4 (parspace explorer) would be served by
% setting a few more options. These settings can be ignored for other
% optimisation routines, but skip_sg must be defined before running
% automatic_runs.
opt_optim.ps_saved = 0; % use saved set for parameter-space explorer (1) or not (0);
opt_optim.ps_plots = 0; % when set to 1, makes intermediate plots of parameter space to monitor progress
opt_optim.ps_profs = 0; % when set to 1, makes profiles and additional sampling for parameter-space explorer
opt_optim.ps_rough = 1; % set to 1 for rough settings of parameter-space explorer, 0 for settings as in openGUTS
skip_sg            = 0; % set to 1 to skip startgrid completely (use ranges in <par> structure)

% The switch fit_tox is used to select whether to fit the control treatment
% (0), the toxicity treatments (1; the control is shown as well), or both
% together (2). A good strategy is to fit the basic parameters to the
% control data first (fit_tox=0), copy the best values into the parameter
% matrix above, and then fit the toxicity parameter to the complete data
% set (fit_tox=1). The code below automatically does that.

% ===== FITTING CONTROLS ==================================================
% Simplex fitting works fine for control data. Note that hb is fitted along
% with the other parameters.

opt_optim.type = 1; % optimisation method: 1) default simplex, 4) parspace explorer
fit_tox = [0]; % fit control parameters in data set
par_out = automatic_runs_basic(fit_tox,par,ind_tox,skip_sg,opt_optim,opt_plot);
% script to run the calculations and plot, automatically
par = copy_par(par,par_out,1); % copy fitted parameters into par, and keep fit mark in par

% ===== FITTING TOX DATA ==================================================
% Here, we can also use the parameter-space explorer for fitting the
% treatments. For the explorer, note that, with fit_tox=1, the tox
% parameters in par are replaced (when fitted) with estimates based on the
% data set using startgrid_debtox (called in automatic_runs). This is still
% quite experimental, so you may need to restart with manually-adapted
% ranges! Furthermore, this is really slow ... (the parallel toolbox really
% helps here!)
% 
% For this demo, I use the regular simplex fitting. The parameter-space
% explorer is demonstrated for the ERA-special package

opt_optim.type = 1; % optimisation method 1) simplex, 4 parameter-space explorer
fit_tox = [1]; % fit tox parameters in data set
par_out = automatic_runs_basic(fit_tox,par,ind_tox,skip_sg,opt_optim,opt_plot);
disp_settings_debtox2019 % display some information on the settings on screen

% dedicated TKTD plots; these plots are more readable, especially when plotting CIs
opt_tktd.repls = 0; % plot individual replicates (1) or means (0)
if opt_optim.type == 4    % when the parspace explorer was used, we can plot with CIs ...
    opt_conf.type    = 3; % make intervals from parspace explorer
    opt_conf.lim_set = 2; % use limited set of n_lim points (1) or outer hull (2) to create CIs
    plot_tktd(par_out,opt_tktd,opt_conf);
else
    plot_tktd(par_out,opt_tktd,[]); % leave the options opt_conf empty to suppress all CIs for these plots
    % Note: plotting CIs requires a sample as saved by the various
    % methods available in BYOM (see examples directory).
end

return % stops the analysis after optimisation/plotting; run the additional analyses below manually (or remove this line)

%% Profiling the likelihood
% By profiling you make robust confidence intervals for one or more of your
% parameters. Use the names of the parameters as they occurs in your
% parameter structure _par_ above. This can be a single string (e.g.,
% 'kd'), a cell array of strings (e.g., {'kd','mw'}), or 'all' to profile
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
% calc_proflik(par_out,'kd',opt_prof,opt_optim);  % calculate a profile
% % Enter single parameter names, a cell array of names, or 'all' to profile
% % all fitted parameters.

%% Calculate ECx for all endpoints
% The function calc_ecx will calculate ECx values, based on the parameter
% values (and a sample from parameter space to make CIs). The ECx is
% calculated assuming constant exposure (which is in its definition) at the
% time points requested.

opt_conf.type    = 0; % make intervals from 1) slice sampler, 2)likelihood region, 3) parspace explorer
opt_conf.lim_set = 1; % use limited set of n_lim points (1) or outer hull (2, not for Bayes) to create CIs
opt_ecx.statsup  = []; % states to suppress from the calculations (e.g., glo.locS)

% % UNCOMMENT FOLLOWING LINE(S) TO CALCULATE
% calc_ecx(par_out,[10 14 21 30],opt_ecx,opt_conf);
% % Note that first input can be par_out, but if left empty it is read from
% % the saved file (if it exists for the specified opt_conf.type). Second
% % input is a vector with time points on which ECx is calculated (may extend
% % beyond the time vector of the data set).

%% Population growth rate
% The function calc_pop allows for a calculation of the intrinsic rate of
% population increase (Euler-Lotka equation). The function calc_pop needs
% checking ... so use with care! By default, the calculation is based on
% continuous reproduction (like the fit).

% Population calculations can be performed with CIs
opt_conf.type    = 0; % make intervals from 1) slice sampler, 2)likelihood region, 3) parspace explorer
opt_conf.lim_set = 0; % for lik-region sample: use limited set of n_lim points (1) or outer hull (2) to create CIs

Tpop = linspace(0,50,100); % time vector for the population calculations
Cpop = linspace(0,1,100); % concentration range for population calculations
% Cpop = -1; % use only concentrations as given in X0mat
Th = 3; % time from fresh egg to t=0 in time vector for experimental test (e.g., hatching time)
% Here, use 3 days, as this value was earlier subtracted from the time
% vector in the data set. So, reproduction now represents egg formation,
% and we need to account for hatching time in the population calculation.

% % UNCOMMENT FOLLOWING LINE(S) TO CALCULATE
% calc_pop(par_out,X0mat,Tpop,Cpop,Th,opt_conf,opt_pop)

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