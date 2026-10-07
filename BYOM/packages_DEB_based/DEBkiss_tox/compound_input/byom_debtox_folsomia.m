%% BYOM, simple model: byom_debtox_folsomia.m
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
% *The model:* DEBtox model for toxicants, based on DEBkiss. Internally, it
% is formulated in primary parameters while the user interacts with more
% friendly compound ones. The model includes flexible modules for
% toxicokinetics/damage dynamics and toxic effects. The DEBkiss e-book (see
% <http://www.debtox.info/book_debkiss.html>) provides a partial
% description of the model; the publication of Jager in Ecological
% Modelling contains the full details:
% <https://doi.org/10.1016/j.ecolmodel.2019.108904>.
%
% *This script:* The springtail _Folsomia candida_ exposed to chlorpyrifos
% in food, data set also used as case study in Jager et al (2007; there
% fitted with a receptor-kinetics module)
% <http://dx.doi.org/10.1016/j.envpol.2006.04.028> and in Jage (2020) with
% a simplified model <https://doi.org/10.1016/j.ecolmodel.2019.108904>.
% Simultaneous fit on growth, reproduction and survival.
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
DATA{1} = [0]; % there are never data for this state

% wet body weight (mg) on each observation time (d), concentrations
% in mg/kg food. Removed the observations beyond t=51 days here.
DATA{2} = [0.5	0	0	0.93	2	4.31	9.28	20
    0   NaN NaN NaN NaN NaN NaN NaN
    11	0.034435714	0.031236364	0.031653846	0.030528571	0.032414286	0.03024	0.032857143
    16	0.071646667	0.066173333	0.077033333	0.063886667	0.06492	0.06834	0.06484
    23	0.141813333	0.120026667	0.12588	0.127233333	0.115426667	0.12238	0.119661538
    30	0.197253333	0.175433333	0.17492	0.19766	0.165713333	0.15356	0.166169231
    37	0.220106667	0.207553333	0.196426667	0.203186667	0.18576	0.192426667	0.191833333
    44	0.233993333	0.217633333	0.213446667	0.22014	0.201006667	0.21082	0.206666667
    51	0.260473333	0.223593333	0.228893333	0.23876	0.226526667	0.216653333	0.215209091
%     58	0.26656	0.235986667	0.236578571	0.249026667	0.215546667	0.22774	0.223554545
%     65	0.276306667	0.265233333	0.253357143	0.265046667	0.25124	0.248173333	0.249081818
%     79	0.276706667	0.27322	0.270514286	0.28732	0.261893333	0.261938462	0.24556
%     93	0.287323077	0.280353846	0.2852	0.309384615	0.284664286	0.278958333	0.2519
%     107	0.273115385	0.26517	0.296146154	0.306445455	0.293938462	0.2788125	0.2693
%     121	0.29285	0.28484	0.309818182	0.282125	0.290608333	0.2887	0.2779
%     138	0.306171429	0.29865	0.32534	0.296733333	0.29668	0.30742	0.286257143
    ];

DATA{2}(2:end,2:end) = (DATA{2}(2:end,2:end)).^(1/3); % mg to mm volumetric length (assume density 1 mg_wet/mm3)

% weights: number of individuals measured
W{2} = [0 0 0 0 0 0 0
    14	11	13	14	14	15	14
    15	15	15	15	15	15	15
    15	15	15	15	15	15	13
    15	15	15	15	15	15	13
    15	15	15	15	15	15	12
    15	15	15	15	15	15	12
    15	15	15	15	15	15	11
%     15	15	14	15	15	15	11
%     15	15	14	15	15	15	11
%     15	15	14	15	15	13	10
%     13	13	14	13	14	12	7
%     13	10	13	11	13	8	7
%     12	5	11	8	12	4	6
%     7	4	10	6	10	5	7
    ]; 

% cumulative reproduction (nr. eggs per individual) on each observation
% time (d), concentrations in mg/kg food
DATA{3} = [0.5	0	0	0.93	2	4.31	9.28	20
    0	NaN NaN NaN NaN NaN NaN NaN
    8	0	0	0	0	0	0	0
    10	0	0	0	0	0	0	0
    12	0	0	0	0	0	0	0
    15	0	0	0	0	0	0	0
    17	11.9	10.1	6.5	2.6	0	0	0
    19	29.8	20.6	18.4	21.3	14.3	2.875	0
    22	29.8	20.6	23.6	21.3	19	6	0
    24	55.8	61	46.2	53.8	25.7	29.57142857	1.714285714
    26	116.4	126.7	114.5	107	94.4	34.19642857	1.714285714
    29	119.4	129.6	128.5	117.7	107.1	51.82142857	1.714285714
    31	119.4	134.6	142.1	117.7	108.4333333	55.67857143	12.04761905
    34	235.3	241.3	236	208.3	177.5444444	90.39285714	12.04761905
    36	238.4	241.3	252.7	231	228.1	116.6785714	12.04761905
    38	241.4	243	252.7	231	230.4333333	116.6785714	12.04761905
    40	248.6	243	252.7	231	230.4333333	132.1071429	12.04761905
    43	312.4	306.3333333	301.5	276.8888889	292.9888889	152.3571429	12.04761905
    45	333.775	334.7777778	327.7	337.031746	329.2388889	188.6904762	12.04761905];

% weights: number of individuals from which reproduction was determined
W{3} = [10	10	10	10	10	10	10
    10	10	10	10	10	10	10
    10	10	10	10	10	10	10
    10	10	10	10	10	10	7
    10	10	10	10	10	8	8
    10	10	10	10	10	8	7
    10	10	10	10	10	8	8
    10	10	10	10	10	8	6
    10	10	10	10	10	7	7
    10	10	10	10	10	8	7
    10	10	10	10	10	8	7
    10	10	10	10	9	7	6
    10	10	10	10	9	7	6
    10	10	10	10	9	7	3
    10	10	10	10	9	6	2
    10	10	10	9	9	7	2
    10	9	10	9	9	8	2
    8	9	10	7	8	6	2]; 

% survivors on each observation time (d), concentrations in mg/kg food
DATA{4} = [-1	0	0	0.93	2	4.31	9.28	20
    0	30	30	30	30	30	30	30
    8	30	30	30	30	30	30	30
    10	30	30	30	30	30	30	30
    12	30	30	30	30	30	30	23
    15	30	30	30	30	30	24	19
    17	30	30	30	30	30	22	14
    19	30	30	30	30	30	22	14
    22	30	30	30	30	30	21	11
    24	30	30	30	30	30	20	11
    26	30	29	30	30	30	20	9
    29	30	28	30	30	30	18	8
    31	30	28	30	30	30	17	6
    34	30	28	30	30	30	17	6
    36	30	27	30	30	30	17	3
    38	30	27	30	30	29	17	2
    40	30	26	30	29	29	17	2
    43	30	25	29	28	29	17	2
    45	28	25	29	26	28	17	2];

% if weight factors are not specified, ones are assumed in start_calc.m

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios.

X0mat = [0 0.93 2 4.31 9.28 20]; % the scenarios (here nominal concentrations) 
X0mat(2,:) = 0; % initial values state 1 (scaled damage)
X0mat(3,:) = 0; % initial values state 2 (body length, initial value overwritten by L0M)
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
% conversion factors for Folsomia candida, and yield factors (keep fixed)
glo.delM  = 1;    % shape corrector (-) (data here are volumetric length)
glo.dV    = 0.28; % dry weight density (mg/mm3)
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
% glo.moa = [0 0 1 1 0]; % costs for growth and reproduction
glo.moa = [0 0 0 1 0]; % costs for reproduction
% glo.moa = [0 0 0 0 1]; % hazards for reproduction

% select which feedbacks to use on damage dynamics as set of switches:
% [surface:volume on uptake, surface:volume on elimination, growth dilution, losses with reproduction] 
% glo.feedb = [1 1 1 1]; % all feedbacks 
% glo.feedb = [1 1 1 0]; % classic DEBtox (no losses with repro)
% glo.feedb = [0 0 1 0]; % damage that is diluted by growth
glo.feedb = [0 0 0 0]; % damage that is not diluted by growth

% syntax: par.name = [startvalue fit(0/1) minval maxval optional:log/normal scale (0/1)];
par.L0M  = [0.122   0 0 1e6 1]; % body length at start experiment (mm)
par.LpM  = [0.4     1 0 1e6 1]; % body length at puberty (mm)
par.LmM  = [0.7     1 0 1e6 1]; % maximum body length (mm)
par.rB   = [0.05    1 0 1e6 1]; % von Bertalanffy growth rate constant (1/d)
par.Rm   = [22      1 0 1e6 1]; % maximum reproduction rate (#/d)
par.WB0  = [0.00041 0 0 1e6 1]; % initial dry weight of egg (mg)
par.f    = [1       0 0   2 1]; % scaled functional response
par.hb   = [0.002   1 0 1e6 1]; % background hazard rate (d-1)
par.a    = [1     0 0.1  10   0]; % coefficient for Weibull backgound hazard (-)
% Note: it does not matter whether length measures are entered as actual
% length or as volumetric length (as long as the same measure is used consistently).

% extra parameters for special situations
par.LfM  = [0 0 0 1e6 1]; % actual body length at half-saturation of feeding (zero to ignore)  (mm)
par.LjM  = [0 0 0 1e6 1]; % actual body length at end acceleration (mm)
par.Tlag = [0 0 0 1e6 1]; % lag time for start development

ind_tox = length(fieldnames(par))+1; % index where tox parameters start
% the parameters below this line are all treated as toxicity parameters!
par.kd   = [0.23   1 0  20 1]; % dominant rate constant (d-1)
par.zb   = [8.6    1 0 1e6 1]; % effect threshold energy budget ([C])
par.bb   = [1.9    1 0 1e6 1]; % effect strength energy-budget effect (1/[C])
par.zs   = [6.7    1 0 1e6 1]; % effect threshold survival ([C])
par.bs   = [0.0047 1 0 1e6 1]; % effect strength survival (1/([C] d))

% you can copy-paste fitted basic parameters below this line


%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% constructed, based on the data set.

% specify the y-axis labels for each state variable
glo.ylab{1} = 'scaled damage mg/kg food)';
glo.ylab{2} = 'volumetric body length (mm)';
if isfield(glo,'Tbp') && glo.Tbp > 0
    glo.ylab{3} = ['cumul. repro. (shift ',num2str(glo.Tbp),'d)'];
else
    glo.ylab{3} = 'cumul. repro. (no shift)';
end
glo.ylab{4} = 'survival fraction (-)';

% specify the x-axis label (same for all states)
glo.xlab    = 'time (days)';
glo.leglab1 = 'conc. '; % legend label before the 'scenario' number
glo.leglab2 = 'mg/kg food'; % legend label after the 'scenario' number

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

opt_optim.type   = 1; % optimisation method: 1) default simplex, 4) parspace explorer
opt_optim.fit    = 1; % fit the parameters (1), or don't (0)
opt_optim.it     = 1; % show iterations of the optimisation (1, default) or not (0)
opt_plot.bw      = 1; % if set to 1, plots in black and white with different plot symbols
opt_plot.annot   = 2; % annotations in multiplot for fits: 1) box with parameter estimates 2) single legend
opt_plot.statsup = [1]; % vector with states to suppress in plotting fits
opt_plot.limax   = 1; % if set to 1, limit axes to the data set for each stage

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