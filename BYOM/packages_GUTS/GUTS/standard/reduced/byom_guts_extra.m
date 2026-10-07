%% BYOM, byom_guts_extra.m
%
% * Author: Tjalling Jager
% * Date: November 2021
% * Web support: <http://www.debtox.info/byom.html>
% * Back to index <walkthrough_guts.html>
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
% reduced model (TK and damage dynamics lumped): SD, IT and mixed (or GUTS
% proper). For faster calculations, scaled damage is calculated using the
% analytical solution in <simplefun.html simplefun.m>. Most types of
% time-varying exposure can use (partial) analytical solutions as well (not
% for continuous splines). Calculation of the death mechanisms (SD, IT and
% mixed) is done in <call_deri.html call_deri.m>.
%
% The equation for scaled damage (referenced to water) is used both by SD
% and IT:
%
% $$ \frac{dD_w}{dt} = k_d \left(C_w - D_w \right) \quad \textrm{scaled damage (referenced to external concentration)} $$
%
% The SD model uses the hazard rate, which is calculated as:
%
% $$ h_z = b_w \max (0,D_w-z_w) + h_b \quad \textrm{with $z_w=m_w$, hazard rate}$$
%
% The hazard rate is integrated over time to yield the survival
% probability:
%
% $$ S_z = \textrm{exp} \left(-\int_0^t{h_z(\tau)}d\tau \right) $$
%
% The IT model first finds the maximum of the scaled damage over time:
%
% $$ D_{wm} = \max_{0<\tau<t} D_w(\tau) $$
%
% And then uses the cumulative distribution of the thresholds to find the
% survival probability:
%
% $$ S = 1-F(D_{wm}|m_w,\beta) \quad \textrm{survival probability} $$
%
% *This script:* Long acute toxicity test for guppy (_Poecilia reticulata_)
% exposed to the insecticide dieldrin. Data from: Bedaux and Kooijman
% (1994), <http://dx.doi.org/10.1007/BF00469427>. In this script, many
% options for post-calculations are provided. You can copy the parts of
% code that you need to the other examples (and your own analyses). The
% files in this directory use the analytical solution only!
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

% observed number of survivors, time in days, conc. in ug/L
DATA{1} = [ -1  0  3.2  5.6  10  18  32  56  100
		0     20   20   20  20  20  20  20   20
		1     20   20   20  20  18  18  17    5
		2     20   20   19  17  15   9   6    0
		3     20   20   19  15   9   2   1    0
		4     20   20   19  14   4   1   0    0
		5     20   20   18  12   4   0   0    0
		6     20   19   18   9   3   0   0    0
		7     20   18   18   8   2   0   0    0 ];

% % demonstration of how static-renewal can be used in BYOM    
% Cw = [  -1 3.2  5.6  10  18  32  56  100 % these are now the scenario identifiers
%          0 3.2  5.6  10  18  32  56  100 % initial concentrations at t=0
%          2 3.2  5.6  10  18  32  56  100 % renewals at t=2
%          4 3.2  5.6  10  18  32  56  100 % renewals at t=4
%          6 3.2  5.6  10  18  32  56  100 % renewals at t=6
%          7 0.2  0.2  0.2 0.2 0.2 0.2 0.2]; % last time point and disappearance rate constant for each scenario
% 
% make_scen(3,Cw); % create the global definition for the exposure scenarios in Cw
      
% Note: For survival data, the weights matrix W has a different meaning
% than for continuous data: it can be used to specify the number of animals
% that went missing or were removed during the experiment (enter the number
% of missing/removed animals at the time they were last seen alive in the
% test). By default, the weights matrix for survival data is filled with
% zeros.
    
% scaled damage, can have no observations
DATA{2} = 0;   

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios.

X0mat(1,:) = DATA{1}(1,2:end); % scenarios (concentrations)
X0mat(2,:) = 1;                % initial survival probability
X0mat(3,:) = 0;                % initial scaled damage

%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. 

% global parameters for GUTS purposes
glo.sel  = 1; % select death mechanism: 1) SD 2) IT 3) mixed
glo.locS = 1; % location of survival probability in the state variable list
glo.locD = 2; % location of scaled damage in the state variable list

% start values for SD
% syntax: par.name = [startvalue fit(0/1) minval maxval optional:log/normal scale (0/1)];
par.kd = [0.5  1 1e-3  100 0];   % dominant rate constant, d-1
par.mw = [5    1    0  1e6 1];   % median threshold for survival (ug/L)
par.hb = [0.01 1    0    1 1];   % background hazard rate (1/d)
par.bw = [0.05 1 1e-6  1e6 0];   % killing rate (L/ug/d) (SD only)
par.Fs = [3    1    1  100 1];   % fraction spread of threshold distribution (IT only)

% % start values for IT
% % syntax: par.name = [startvalue fit(0/1) minval maxval optional:log/normal scale (0/1)];
% par.kd    = [0.008 1 1e-3 100 0]; % dominant rate constant (d-1)
% par.mw    = [0.5   1    0 1e6 1]; % median threshold for survival (ug/L)
% par.hb    = [0.001 1    0   1 1]; % background hazard rate (1/d)
% par.bw    = [0.04  1 1e-6 1e6 0]; % killing rate (L/ug/d) (SD and mixed)
% par.Fs    = [3.8   1    1 100 1]; % fraction spread of NEC distribution (-) (IT and mixed)

switch glo.sel % make sure that right parameters are fitted
    case 1 % for SD ...
        par.Fs(2) = 0; % never fit the threshold spread!
    case 2 % for IT ...
        par.bw(2) = 0; % never fit the killing rate!
    case 3 % mixed
        % do nothing: fit all parameters
end

%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% constructed, based on the data set

% specify the y-axis labels for each state variable
glo.ylab{1} = 'survival probability';
glo.ylab{2} = ['scaled damage (',char(181),'g/L)'];
% specify the x-axis label (same for all states)
glo.xlab    = 'time (days)';
glo.leglab1 = 'conc. '; % legend label before the 'scenario' number
glo.leglab2 = [char(181),'g/L']; % legend label after the 'scenario' number

prelim_checks % script to perform some preliminary checks and set things up
% Note: prelim_checks also fills all the options (opt_...) with defauls, so
% modify options after this call, if needed.

%% Calculations and plotting
% Here, the function is called that will do the calculation and the
% plotting. Options for the plotting can be set using opt_plot (see
% prelim_checks.m). Options for the optimisation routine can be set using
% opt_optim. The files in this directory always apply the analytical
% solution for damage in simplefun.
% 
% For the demo, the iterations are turned off (opt_optim.it = 0).

% par = start_vals_guts(par); % experimental start-value finder; use at your own risk
% % Note: start_vals will now overwrite the parameter structure par!

% glo.fastslow = 's'; % tell simplefun that we need slow (s) or fast (f) kinetics, (o) is off (default)
% par.kd(2)    = 0; % don't fit the dominant rate constant anymore
% % slow kinetics implies that mw is now used for mw/kd and bw becomes bw*kd!
% par.mw = [60  1 0 1e6 1]; % median threshold-time for survival (d ug/L)
% % here, I only define a new starting value for mw/kd, assuming we're fitting IT

opt_optim.fit    = 1; % fit the parameters (1), or don't (0)
opt_optim.it     = 1; % show iterations of the optimisation (1, default) or not (0)
opt_plot.bw      = 1; % if set to 1, plots in black and white with different plot symbols
opt_plot.annot   = 2; % annotations in multiplot for fits: 1) box with parameter estimates 2) single legend
% opt_plot.statsup = [2]; % vector with states to suppress in plotting fits

% optimise and plot (fitted parameters in par_out)
par_out = calc_optim(par,opt_optim); % start the optimisation
calc_and_plot(par_out,opt_plot); % calculate model lines and plot them

opt_tktd.obspred = 0; % plot predicted-observed plots (1) or not (0)
plot_tktd(par_out,opt_tktd,[]); % make more detailed standard TKTD plots (without CIs)

% savenm = 'fit_dieldrin'; % file name to save a plot ad hoc
% save_plot(gcf,savenm,[],3) % save active figure as PDF (see save_plot for more options)
% print_par(par_out) % print fitted parameters on-screen in format that can be copied into this script

% Below, several options for post analyses are provided.

return % stop here, you can run the next sections afterwards if you like

%% Local sensitivity analysis
% Local sensitivity analysis of the model. All model parameters are
% increased one-by-one by a small percentage. The sensitivity score is
% by default scaled (dX/X p/dp) or alternatively absolute (dX p/dp).
%
% Options for the sensitivity can be set using opt_sens (see
% prelim_checks.m).
  
% % UNCOMMENT FOLLOWING LINE(S) TO CALCULATE
% calc_localsens(par_out,opt_sens)

%% Profiling the likelihood
% By profiling you make robust confidence intervals for one or more of your
% parameters. Use the names of the parameters as they occurs in your
% parameter structure _par_ above. This can be a single string (e.g.,
% 'kd'), a cell array of strings (e.g., {'kd','mw'}), or 'all' to profile
% all fitted parameters. This example produces a profile for each parameter
% and provides the 95% confidence interval (on screen and indicated by the
% horizontal broken line in the plot).
%
% *Note: for more post-calculations, see <byom_bioconc_extra.html
% byom_bioconc_extra.m>*
%
% Options for profiling can be set using opt_prof (see prelim_checks.m).
% 
% For this demo, no sub-optimisations are used. However, consider this
% options when working on your own data (e.g., set opt_prof.subopt=10). The
% level of detail of the profiling is set to 'detailed' for this demo to
% catch the local optimum (it is just within the 95% CI, which is thus a
% broken set)

opt_prof.detail   = 1; % detailed (1) or a coarse (2) calculation
opt_prof.subopt   = 10; % number of sub-optimisations to perform to increase robustness
opt_prof.brkprof  = 2; % when a better optimum is located, stop (1) or automatically refit (2)

% % UNCOMMENT FOLLOWING LINE(S) TO CALCULATE
% calc_proflik(par_out,'all',opt_prof,opt_optim);  % calculate a profile
% % Enter single parameter names, a cell array of names, or 'all' to profile
% % all fitted parameters.

%% Slice sampler
% The slice sampler can be used for a Bayesian analysis as it provides a
% sample from the posterior distribution. A .mat file is saved which
% contains the sample, to use later for e.g., intervals on model
% predictions. The output includes the Markov chain and marginal
% distributions for each fitted parameter. The sample can be used to put
% confidence intervals on the model lines, as demonstrated below.
%
% Options for the slice sampling can be set using opt_slice (see
% prelim_checks). Options for confidence bounds on model curves and
% calculation of error ellipse can be set using opt_conf (see
% prelim_checks).
%
% For this data set, it is difficult to obtain a proper posterior (one that
% can be integrated so that sample can be normalised to a prob. density
% function). For SD, _hb_ runs towards zero, and for IT, _kd_ runs to zero.
% Some careful selection of priors will be needed to make the posterior
% proper. Therefore, this analysis is skipped here, and only the
% likelihood-region method is shown.

% % UNCOMMENT FOLLOWING LINE(S) TO CALCULATE 
% opt_slice.thin     = 25; % thinning of the sample (keep one in every 'thin' samples)
% opt_slice.burn     = 200; % number of burn-in samples (0 is no burn in)
% opt_slice.alllog   = 0; % set to 1 to put all parameters on log-scale before taking the sample
% 
% calc_slice(par_out,1000,opt_slice); % second argument number of samples (-1 to re-use saved sample from previous runs)

%% Likelihood region (shooting)
% Another way to make intervals on model predictions is to use a sample of
% parameter sets taken from the joint likelihood-based conf. region. This
% is done by the function calc_likregion.m. It first does profiling of all
% fitted parameters to find the edges of the region. Then, Latin-Hypercube
% shooting, keeping only those parameter combinations that are not rejected
% at the 95% level in a lik.-rat. test. The inner rim will be used for CIs
% on forward predictions.
% 
% In general, it is a good idea to take 10 sub-optimisations for
% robustness. In this case, skipping sub-optimisation would mean that the
% profile for _hb_ is not fully accurate. Options for the likelihood region
% can be set using opt_likreg (see prelim_checks.m).

opt_prof.detail  = 1;  % detailed (1) or a coarse (2) calculation
opt_prof.subopt  = 10; % number of sub-optimisations to perform to increase robustness
opt_prof.brkprof = 2; % when a better optimum is located, stop (1) or automatically refit (2)
opt_likreg.burst = 2000; % number of random sets from parameter space evaluated every iteration
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
% Plotting with sampling error is restricted to Bayesian analyses for the
% moment (for the frequentist context, it needs more work), and only with
% plot_guts (not demonstrated). Options for confidence bounds on model
% curves can be set using opt_conf (see prelim_checks).
% 
% Use opt_conf.type to tell calc_conf which sample to use: 
% -1) Skips CIs (zero does the same, and an empty opt_conf as well).
% 1) Bayesian MCMC sample (default); CI on predictions are 95% ranges on 
% the model curves from the sample 
% 2) parameter sets from a joint likelihood region using the shooting 
% method (limited sets can be used), which will yield (asymptotically) 95% 
% CIs on predictions
% 3) as option 2, but using the parameter-space explorer
%
% The plot_tktd function makes multiplots for the survival data, which are
% more readable when plotting with various intervals. It plots survival and
% damage versus time, with each treatment in a separate panel. Furthermore,
% a predicted-observed plot is made. CIs on observed survival probabilities
% are the Wilson score intervals. Note that plot_tktd calls calc_conf
% itself, so there is no need to run it here first to obtain out_conf. The
% options opt_conf are also used by plot_tktd. The plot_tktd function is
% meant to replace plot_guts (which is still available in the engine, if
% needed).

opt_conf.type    = 2; % make intervals from 1) slice sampler, 2)likelihood region shooting, 3) parspace explorer
opt_conf.lim_set = 2; % for type 2/3: use limited set of n_lim points (1) or outer hull (2) to create CIs
opt_tktd.obspred = 1; % plot predicted-observed plots (1) or not (0)

% % Make standard plot with CIs
% out_conf = calc_conf(par_out,opt_conf,opt_plot.repls); % calculate confidence intervals on model curves
% % Note: extra argument opt_plot.repls is needed to obtain correct sampling error (so when opt_conf.samerr = 1, and Bayes only).
% calc_and_plot(par_out,opt_plot,out_conf); % call the plotting routine again to make fits with CIs (not so readable)

% Calculate and plot specific plots for TKTD models, more readable with
% CIs. This function replaces plot_guts.
plot_tktd(par_out,opt_tktd,opt_conf); 

%% Calculate LCx versus time
% Here, the LCx (by default the LC50) is calculated at several time points.
% When sufficient points are specified, a smooth line for LCx versus time
% will be produced. LCx values are also printed on screen. If a sample from
% parameter space is available (e.g., from the slice sampler or the
% likelihood region), it can be used to calculate confidence bounds. Note
% that opt_conf.type=-1 skips CIs.
% 
% Here, a faster method is used with calc_lcx_lim, which only works with
% the standard models as provided in this package (and in this case only
% the reduced model). As soon as you modify the model (in call_deri,
% derivatives or simplefun) it will likely produce erroneous results. The
% function calc_ecx is always applicable, but slower.
%
% Options for LCx (with confidence bounds) can be set using opt_lcx_lim
% (see prelim_checks). Note that the calculations assume initial values for
% each state as defined in the first column of X0mat! Note that
% opt_conf.type=-1 skips CIs.

opt_conf.type    = 2; % make intervals from 1) slice sampler, 2)likelihood region shooting, 3) parspace explorer
opt_conf.lim_set = 2; % for lik-region sample: use limited set of n_lim points (1) or outer hull (2) to create CIs
opt_lcx_lim.Feff = 0.50; % effect level (>0 en <1), x/100 in LCx/LPx

Tend = [1:8]; % times at which to calculate LCx, relative to control
calc_lcx_lim_guts_red(par_out,Tend,opt_lcx_lim,opt_conf); % calculates LCx values, CI requires that there is a mat file with sample

% % This is the slower general method
% opt_ecx.Feff      = [0.10 0.50]; % effect levels (>0 en <1), x/100 in ECx
% opt_ecx.notitle   = 1; % set to 1 to suppress titles above ECx plots
% Tend = [1:8]; % times at which to calculate LCx, relative to control
% calc_ecx(par_out,Tend,opt_ecx,opt_conf); % general method for ECx values

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
%
% This is the fast method, which is closely linked to this GUTS model
% (reduced model, standard directory). Don't use if you modify the GUTS
% model in simplefun.m, call_deri.m or derivatives.m, unless you know what
% you are doing! The function calc_epx can also calculate LPx values. That
% function works in all cases, but is much slower.

opt_conf.type     = 2; % make intervals from 1) slice sampler, 2)likelihood region shooting, 3) parspace explorer
opt_conf.lim_set  = 2; % for lik-region sample: use limited set of n_lim points (1) or outer hull (2) to create CIs
opt_lcx_lim.Feff  = 0.10; % effect level (>0 en <1), x/100 in LCx/LPx

Tend  = []; % time at which to calculate LPx (empty is end of profile)
LPx_mon = calc_lpx_lim_guts_red(par_out,Tend,'profile_monit.txt',opt_lcx_lim,opt_conf);

% % This is the slower general method
% opt_ecx.Feff      = [0.10]; % effect levels (>0 en <1), x/100 in ECx
% opt_ecx.notitle   = 1; % set to 1 to suppress titles above ECx plots
% calc_epx(par_out,'profile_monit.txt',[],opt_ecx,opt_conf,opt_tktd); % general method for EPx values

%% Other files: simplefun
% To archive analyses, publishing them with Matlab is convenient. To keep
% track of what was done, the file simplefun.m can be included in the
% published result.
% 
% <include>simplefun.m</include>

%% Other files: call_deri
% To archive analyses, publishing them with Matlab is convenient. To keep
% track of what was done, the file call_deri.m can be included in the
% published result.
%
% <include>call_deri.m</include>