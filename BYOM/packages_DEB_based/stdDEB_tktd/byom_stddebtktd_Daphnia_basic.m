%% BYOM, byom_stddebtktd_Daphnia_basic.m
%
% * Author: Tjalling Jager
% * Date: June 2022
% * Web support: <http://www.debtox.info/byom.html>
%
% BYOM is a General framework for simulating model systems in terms of
% ordinary differential equations (ODEs). The model itself needs to be
% specified in <derivatives.html derivatives.m>, and <call_deri.html
% call_deri.m> may need to be modified to the particular problem as well.
% The files in the engine directory are needed for fitting and plotting.
% Results are shown on screen but also saved to a log file (results.out).
%
% *The model:* standard DEB model in terms of powers (energy fluxes in J/d)
% and body length. In this package, toxicant stress is included. The model
% equations are those from the 'Tromso coffee mug' ;-). Only change is that
% the ODE for structure is phrased in volumetric length rather than volume.
% The TKTD model is lifted from DEBtox2019. The 2023 publication of Jager
% et al in Ecological Modelling contains the full details:
% <https://doi.org/10.1016/j.ecolmodel.2022.110187>.
%
% *This script:* The water flea _Daphnia magna_ exposed to fluoranthene;
% data set from Jager et al (2010),
% <http://dx.doi.org/10.1007/s10646-009-0417-z>. Published as case study in
% Jager & Zimmer (2012). Simultaneous fit on growth, reproduction and
% survival. Here the mean data from 5-10 females is used. This means that
% the special functions for censoring the reproduction data cannot be used.
% This is a simple version, using automatic_runs_basic.
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

% NOTE: there are two controls in the data set, but in this analysis, they
% will have to be lumped. Check byom_stddebtktd_Daphnia_ERA for an example
% where they are separated, and individuals are followed instead of mean
% responses.

% body length (mm) on each observation time (d), concentrations in uM
DATA{3} = [1 0 0 0.213	0.426	0.853
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

DATA{3}(2:end,2:end) = DATA{3}(2:end,2:end) / 10; % from mm to cm

W{3} = 5 * ones(size(DATA{2})-1); % weights: 5 animals in each treatment

% cumulative reproduction (nr. offspring) (mm) on each observation time (d), concentrations in uM
DATA{4} = [1 0 0 	0.213	0.426	0.853
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
W{4} = [10	10	10	10	10
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

% survivors on each observation time (d), concentrations in uM
DATA{6} = [-1	0	0	0.213	0.426	0.853
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

save([glo.basenm,'_DATA'],'DATA','W') % save MAT file with data set 
% Saving the data set is handy to allow for a simple reconstruction of the
% calibrations, without needing to define the data again, in the same way.
% Note that if you use DATAx and Wx, you need to save them as well!

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios. Note that for stdDEB, the initial values cannot
% be randomly chosen. We need to simulate the individual up to the point
% where we start the model analysis (e.g., at birth). Therefore, most of
% the values in X0mat will not be used in the analysis.

X0mat = [0 0.213 0.426 0.853]; % the scenarios (here: exposure concentrations)
X0mat(2,:) = 0;   % initial reserve energy (J) NOT USED
X0mat(3,:) = 0;   % initial maturity level (J) NOT USED
X0mat(4,:) = 0;   % initial volumetric length (cm) NOT USED
X0mat(5,:) = 0;   % initial cumul. repro NOT USED
X0mat(6,:) = 0;   % initial values scaled damage
X0mat(7,:) = 1;   % initial values survival probability

% Put the position of the various states in globals, to make sure correct
% one is selected for extra things (e.g., to prevent shrinking in
% call_deri, and for using plot_tktd).
glo.locL = 3; % location of body size in the state variable list
glo.locR = 4; % location of cumulative reproduction in the state variable list
glo.locD = 5; % location of damage in the state variable list 
glo.locS = 6; % location of survival probability in the state variable list

%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. 
% Global parameters as part of the structure glo

% Global settings for the analysis; several of these values will be
% species-specific. These are not taken from AmP.
glo.yP      = 0.64;  % product of yVE and yEV (-)
glo.len     = 1;     % switch to fit physical length (0=wwt, 1=phys. length, 2=phys. length, no shrinking) (used in call_deri.m)
glo.Tref    = 20 + 273.15; % translate reference temperature to temperature in Kelvin
glo.T       = 20 + 273.15; % translate actual temperature to temperature in Kelvin
glo.Tbp     = 3;     % brood pouch delay (d)
glo.FBV     = 0.02;  % dry weight egg as fraction of structural body weight (-) (for losses with repro; approx. for Daphnia magna)
glo.KRV     = 1;     % part. coeff. repro buffer and structure (kg/kg) (for losses with reproduction)
glo.Lwm_ref = 0.4;   % reference max *physical* length for scaling rate constants (cm)
% Note: using a fixed reference length for scaling is helpful for comparing
% parameters between data sets, and absolutely needed when fitting on
% multiple data sets where animals reach a different maximum length.

glo.E0_calc = [1 1]; % strategy for maternal effects rule
% first element, initial values from provided f in treatment (1) or always use f=1 (2)
% second element, for repro, egg costs from f (1), always use f=1 (2), or from actual reserve status of mother (3)

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

% =========================================================================
%      Parameters for Daphnia magna (Bas Kooijman, Andre Gergs. 2019. 
%      AmP Daphnia magna, version 2019/03/16.)
% =========================================================================
% global settings for conversions from AmP
glo.delM = 0.264; % shape corrector (-)
glo.dV   = 0.19;  % dry weight density (g/cm3)
% syntax: par.name = [startvalue fit(0/1) minval maxval optional:log/normal scale (0/1)];
par.TA   = [6400     0 1000  1e6 1]; % Arrhenius temperature (K)
par.spAm = [313.169  0    0  1e6 1]; % max. surface-specific assimilation rate (J/(cm2 d))
par.spM  = [1200     0    0  1e6 1]; % volume-specific somatic maintenance costs (J/(cm3 d))
par.spT  = [0        0    0  1e6 1]; % surface-specific maintenance costs (J/(cm2 d))
par.kJ   = [0.2537   0    0  1e6 1]; % maturity maintenance rate constant (1/d)
par.EG   = [4400     0    0  1e6 1]; % volume-specific costs for growth (J/cm3)
par.EHb  = [0.05464  0 1e-6  1e6 1]; % maturity level at birth (J)
par.EHj  = [0        0    0  1e6 1]; % maturity level at metamorphosis (J)
par.EHp  = [1.09     1 1e-6  1e6 1]; % maturity level at puberty (J)
par.v    = [0.1858   0    0  1e6 1]; % energy conductance (cm/d)
par.kap  = [0.5809   1 0.01 0.99 1]; % allocation fraction to soma (-)
par.kapR = [0.95     0 0.01 0.99 1]; % reproduction efficiency (-)
par.f    = [1        1    0    2 1];   % scaled food density (-)
par.hb   = [0.01     1 1e-3 0.07 0]; % background hazard rate (d-1)
par.a    = [1        0  0.1   10 0]; % coefficient for Weibull backgound hazard (-) (1 is off)
% par.Lw0  = [0.8      0    0  1e6 1]; % starting at a *physical* length L0>Lb (when glo.len=0, use wet weight)

% NOTE: some parameters for the AmP entry need to be refitted to obtain a
% closer correspondence to the controls of this data set.

% =========================================================================
%      TKTD Parameters for this data set
% =========================================================================
ind_tox = length(fieldnames(par))+1; % index where tox parameters start
% the parameters below this line are all treated as toxicity parameters!

% When using the parameter-space explorer, start values and ranges are not
% needed (filled later by startgrid_debtox); only the fit/fix mark is
% relevant here. But make sure that the value in the first column is within
% the bounds (and not zero for log-scale parameters).
par.kd   = [0.08   1 0.01  10 0]; % dominant rate constant (d-1)
par.zb   = [0.10   1 0    1e6 1]; % effect threshold energy budget ([C])
par.bb   = [72     1 1e-6 1e6 0]; % effect strength energy-budget effect (1/[C])
par.zs   = [0.24   1 0    1e6 1]; % effect threshold survival ([C])
par.bs   = [0.64   1 1e-6 1e6 0]; % effect strength survival (1/([C] d))

%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% used, based on the data set

% specify the y-axis labels for each state variable
glo.ylab{1} = 'reserve (J)';
glo.ylab{2} = 'maturity (J)';
glo.ylab{3} = 'body length (cm)';
if isfield(glo,'Tbp') && glo.Tbp > 0
    glo.ylab{4} = ['cumul. repro (shift ',num2str(glo.Tbp),'d)'];
else
    glo.ylab{4} = 'cumul. repro (no shift)';
end
glo.ylab{5} = ['damage (',char(181),'g/L)'];
glo.ylab{6} = 'survival fraction (-)';
% specify the x-axis label (same for all states)
if isfield(par,'L0')
    glo.xlab    = 'time since start (days)';
else
    glo.xlab    = 'time since birth (days)';
end
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
 
% -------------------------------------------------------------------------
% Configurations for the ODE solver
glo.stiff = [0 3]; % ODE solver 0) ode45 (standard), 1) ode113 (moderately stiff), 2) ode15s (stiff)
% Second argument is for default sloppy (0), normally tight (1), tighter
% (2), or very tight (3) tolerances. Use 1 for quick analyses, but check
% with 3 to see if there is a difference! Especially for time-varying
% exposure, there can be large differences between the settings!
glo.break_time = 0; % break time vector up for ODE solver (1) or don't (0)
% Note: breaking the time vector is a good idea when the exposure scenario
% contains discontinuities. Don't use for continuous splines (type 1) as it
% will be much slower. For FOCUS scenarios (high time resolution), breaking
% up is not efficient and does not appear to be necessary.
% -------------------------------------------------------------------------

opt_optim.fit    = 1; % fit the parameters (1), or don't (0)
opt_plot.bw      = 0; % if set to 1, plots in black and white with different plot symbols
opt_plot.annot   = 2; % annotations in multiplot for fits: 1) box with parameter estimates 2) single legend
opt_plot.repls   = 0; % set to 1 to plot replicates, 0 to plot mean responses

% The switch fit_tox is used to select whether to fit the control treatment
% (0), the toxicity treatments (1; the control is shown as well), or both
% together (2). A good strategy is to fit the basic parameters to the
% control data first (fit_tox=0), copy the best values into the parameter
% matrix above, and then fit the toxicity parameter to the complete data
% set (fit_tox=1). The code below automatically does that.

skip_sg = 0; % set to 1 to skip startgrid completely (use ranges in <par> structure)
% This option is for the parameter-space explorer only.

% ===== FITTING CONTROLS ==================================================
% Simplex fitting works fine for control data. Note that hb is fitted along
% with the other parameters here.
opt_optim.type   = 1; % optimisation method: 1) default simplex, 4) parspace explorer
glo.stdDEB_start = {}; % make sure it is empty and not defined!
fit_tox = [0]; % fit control parameters in data set
par_out = automatic_runs_basic(fit_tox,par,ind_tox,skip_sg,opt_optim,opt_plot);
% script to run the calculations and plot, automatically
par = copy_par(par,par_out,1); % copy fitted parameters into par, and keep fit mark in par

% Next, we need to clear the zero-variate data points, so they are not used
% again for fitting the treatments.
glo.zvd = []; % this should be enough, since it makes the output zvd of call_deri empty as well

% ===== FITTING TOX DATA ==================================================
% Here, we can also use the parameter-space explorer for fitting the
% treatments. Note that, with fit_tox=1, the tox parameters in par are
% replaced (when fitted) with estimates based on the data set using
% startgrid_debtox (called in automatic_runs_basic). This is still quite
% experimental, so you may need to restart with manually-adapted ranges!
% Furthermore, the explorer is really slow ... (the parallel toolbox really
% helps here!). The example below is set up for simplex optimisation, which
% is must faster but not robust unless you have good starting values.

opt_optim.type     = 1; % optimisation method 1) simplex, 4 parameter-space explorer
opt_optim.ps_saved = 0; % use saved set for parameter-space explorer (1) or not (0);
opt_optim.ps_plots = 0; % when set to 1, makes intermediate plots of parameter space to monitor progress
opt_optim.ps_profs = 0; % when set to 1, makes profiles and additional sampling for parameter-space explorer
opt_optim.ps_rough = 1; % set to 1 for rough settings of parameter-space explorer, 0 for settings as in openGUTS

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
