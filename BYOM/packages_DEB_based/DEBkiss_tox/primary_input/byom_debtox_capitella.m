%% BYOM, byom_debkiss_capitella.m
%
% *Table of contents*

%% About
% * Author: Tjalling Jager
% * Date: September 2022
% * Web support: <http://www.debtox.info/byom.html>
%
% BYOM is a General framework for simulating model systems in terms of
% ordinary differential equations (ODEs). The model itself needs to be
% specified in <derivatives.html derivatives.m>, and <call_deri.html
% call_deri.m> may need to be modified to the particular problem as well.
% The files in the engine directory are needed for fitting and plotting.
% Results are shown on screen but also saved to a log file (results.out).
%
% *The model:* DEBtox model for toxicants, based on DEBkiss. This is the
% model formulated in primary parameters. The model includes flexible
% modules for toxicokinetics/damage dynamics and toxic effects. The DEBkiss
% e-book (see <http://www.debtox.info/book_debkiss.html>) provides a
% partial description of the model; the publication of Jager in Ecological
% Modelling contains the full details:
% <https://doi.org/10.1016/j.ecolmodel.2019.108904>.
%
% *This script:* The marine polychaete worm _Capitella teleta_ exposed to
% nonylphenol in sediment. The data set is from Hansen et al (Ecol Appl
% 9:482-295) and has been analysed with standard DEB before by Jager &
% Selck (J Sea Res 66:456-462). Note: here, initial body length is
% calculated from the fixed egg volume.
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
        
% body volume (mm3) on each observation time (d), concentrations in mg/kg sed.
DATA{2} = [1 0	14	52	174
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

% convert volume to volumetric length (mm)
DATA{2}(2:end,2:end) = DATA{2}(2:end,2:end).^(1/3); 

% cumulative reproduction (eggs) on each observation time (d), concentrations in mg/kg sed.
DATA{3} = [1 0	14	52	174
    14	0	0	0	0
    21	0	0	0	0
    25	0	0	0	0
    32	0	0	0	0
    39	493     604     343     0
    46	1007	1234	775     212
    53	1461	1768	1145	427
    60	1936	2284	1599	655
    66	2400	2914	2065	845
    72	2868	3441	2464	1033
    78	3256	3901	2876	1199];

% survival probability (no data)
DATA{4} = 0;

% if weight factors are not specified, ones are assumed in start_calc.m

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios.

X0mat = DATA{2}(1,2:end); % the scenarios (here nominal concentrations) 
X0mat(2,:) = 0; % initial values state 1 (scaled damage)
X0mat(3,:) = 0; % initial values state 2 (body weight, initial value overwritten by calculation from WB0)
X0mat(4,:) = 0; % initial values state 3 (cumulative reproduction)
X0mat(5,:) = 1; % initial values state 4 (survival probability)
         
% X0mat(:,[1]) = []; % remove the control (optional for clearer graphs of the tox data)

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
% conversion factors for Daphnia magna, and yield factors (keep fixed).
glo.delM  = 1;    % shape corrector (-), set to 1 as volume was measured
glo.dV    = 0.2;  % dry weight density (mg/mm3)
glo.yAV   = 0.8;  % yield of assimilates on structure (starvation) (-)
glo.yBA   = 0.95; % yield of egg buffer on assimilates (-)
glo.yVA   = 0.8;  % yield of structure on assimilates (growth) (-)
glo.KRV   = 1;    % part. coeff. repro buffer and structure (kg/kg) (for losses with reproduction)

% settings for the analysis; modify as you like
glo.len   = 1;    % switch to fit length (0=dwt, 1=lenght, 2=length and no shrinking) (used in call_deri.m)
glo.mat   = 1;    % include maturity maint. (0=off, 1=include)

% select mode of action of toxicant as set of switches:
% [assimilation/feeding, maintenance costs (somatic and maturity), growth costs, repro costs] 
% glo.moa = [1 0 0 0 0]; % assimilation/feeding
% glo.moa = [0 1 0 0 0]; % costs for maintenance 
glo.moa = [0 0 1 1 0]; % costs for growth and reproduction
% glo.moa = [0 0 0 1 0]; % costs for reproduction
% glo.moa = [0 0 0 0 1]; % hazards for reproduction

% select which feedbacks to use on damage dynamics as set of switches:
% [surface:volume on uptake, surface:volume on elimination, growth dilution, losses with reproduction] 
% glo.feedb = [1 1 1 1]; % all feedbacks 
glo.feedb = [1 1 1 0]; % classic DEBtox (no losses with repro)
% glo.feedb = [0 0 1 0]; % damage that is diluted by growth
% glo.feedb = [0 0 0 0]; % damage that is not diluted by growth

% syntax: par.name = [startvalue fit(0/1) minval maxval optional:log/normal scale (0/1)];
par.sJAm = [0.151  1 0 1e6 1]; % specific assimilation rate (mg/mm^2/d)
par.sJM  = [0.0514 1 0 1e6 1]; % specific maintenance costs (mg/mm^3/d)
par.kap  = [0.894  1 0   1 1]; % allocation fraction to soma (-)
par.WB0  = [4.1e-3*glo.dV 0 0 1e6 1]; % initial dry weight of egg (see Jager & Selck) (mg)
par.LpM  = [1.84   1 0 1e6 1]; % actual body length at puberty (mm)
par.L0M  = [1      0 0 1e6 1]; % actual body length at start experiment (not used as initial length is calculated from egg weight) (mm)
par.f    = [1      0 0  2 1];  % scaled functional response (-)
par.hb   = [0      0 0 1e6 1]; % background hazard rate (d-1)
par.a    = [1     0 0.1  10   0]; % coefficient for Weibull backgound hazard (-)

% extra parameters for specific cases
par.LfM  = [0      0 0 1e6 1]; % actual body length at half-saturation of feeding (zero to ignore) (mm)
par.LjM  = [0      0 0 1e6 1]; % actual body length at end acceleration (mm)
par.Tlag = [9.00   1 0 30 1];  % lag time for start development (d)
% % alternatively, you can try food limitation to capture the initial growth
% par.LfM  = [0.3    1 0 1e6 1]; % actual body length at half-saturation of feeding (zero to ignore) (mm)
% par.LjM  = [0      0 0 1e6 1]; % actual body length at end acceleration (mm)
% par.Tlag = [0      0 0 30 1];  % lag time for start development (d)
% % or metabolic acceleration
% par.LfM  = [0      0 0 1e6 1]; % actual body length at half-saturation of feeding (zero to ignore) (mm)
% par.LjM  = [0.9    1 0 1e6 1]; % actual body length at end acceleration (mm)
% par.Tlag = [0      0 0 30 1];  % lag time for start development (d)

ind_tox = length(fieldnames(par))+1; % index where tox parameters start
% The parameters below this line are all treated as toxicity parameters!
par.kd   = [0.1   1 0.01  10 0]; % dominant rate constant (d-1)
par.zb   = [4.3   1 0    1e6 1]; % effect threshold energy budget ([C])
par.bb   = [0.007 1 1e-6 1e6 0]; % effect strength energy-budget effects (1/[C])
par.zs   = [500   0 0    1e6 1]; % effect threshold survival ([C])
par.bs   = [0.001 0 1e-6 1e6 0]; % effect strength survival (1/([C] d))

% extra parameter to deal with the hormetic response: f in treatments is a bit higher
par.fh   = [1.05  1 0  2 1]; % hormetic scaled functional response (used in call_deri)

% you can copy-paste fitted basic parameters below this line


%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% used, based on the data set

% specify the y-axis labels for each state variable
glo.ylab{1} = 'scaled damage (mg/kg)';
glo.ylab{2} = 'volumetric body length (mm)';
glo.ylab{3} = 'cumululative reproduction (eggs)';
glo.ylab{4} = 'survival probability';
% specify the x-axis label (same for all states)
glo.xlab    = 'time (days)';
glo.leglab1 = 'conc. '; % legend label before the 'scenario' number
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
glo.break_time = 1; % break time vector up for ODE solver (1) or don't (0)
% Note: breaking the time vector is a good idea when the exposure scenario
% contains discontinuities. Don't use for continuous splines (type 1) as it
% will be much slower. For FOCUS scenarios (high time resolution), breaking
% up is not efficient and does not appear to be necessary.
% 
% Note: here, use glo.break_time=1! Reason is that the lag time is also
% added to the 'events' such that the time vector before the lag is
% run separately. Breaking up the time vector is the best thing to do.
% Alternatively, ode113 has less problems with this lag switch.

opt_optim.type   = 1; % optimisation method: 1) default simplex, 4) parspace explorer
opt_optim.fit    = 1; % fit the parameters (1), or don't (0)
opt_optim.it     = 1; % show iterations of the optimisation (1, default) or not (0)
opt_plot.bw      = 0; % if set to 1, plots in black and white with different plot symbols
opt_plot.annot   = 2; % annotations in multiplot for fits: 1) box with parameter estimates 2) single legend
opt_plot.repls   = 0; % set to 1 to plot replicates, 0 to plot mean responses
opt_plot.statsup = [4]; % vector with states to suppress in plotting fits

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
% calc_proflik(par_out,'all',opt_prof,opt_optim);  % calculate a profile
% % Enter single parameter names, a cell array of names, or 'all' to profile
% % all fitted parameters.

%% Population growth rate
% The function calc_pop allows for a calculation of the intrinsic rate of
% population increase (Euler-Lotka equation). The function calc_pop needs
% checking ... By default, the calculation is based on continuous
% reproduction (like the fit).

% Population calculations can be performed with CIs
opt_conf.type    = 0; % use the values from the slice sampler (1) or likelihood region (2) to make intervals
opt_conf.lim_set = 0; % for lik-region sample: use limited set of n_lim points (1) or outer hull (2) to create CIs

% % UNCOMMENT FOLLOWING LINE(S) TO CALCULATE
% % The extra parameter fh hampers the calculatations for different food
% % levels. So we remove that field here. The population calculations thereby
% % ignore hormesis.
% par_out = rmfield(par_out,'fh');
% 
% Tpop = linspace(0,160,100); % time vector for the population calculations
% Cpop = linspace(0,180,100); % concentration range for population calculations
% % Cpop = -1; % use only concentrations as given in X0mat
% Th = 1; % time from fresh egg to t=0 in time vector for experimental test (e.g., hatching time)
% % Here, it is a guess (based on Ramskov & Forbes, 2008, DOI 10.3354/meps07584).
% 
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