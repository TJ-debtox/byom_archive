%% BYOM, byom_guts_constant_cadmium.m, full GUTS model
%
% * Author: Tjalling Jager 
% * Date: November 2021
% * Web support: <http://www.debtox.info/byom.html>
% * Back to index <walkthrough_guts.html>
%
% BYOM is a General framework for simulating model systems. The files in
% this directory use an analytical solution only, and therefore
% <derivatives.html derivatives.m> will be missing.
%
% *The model:* fitting survival data with the
% <http://www.debtox.info/about_guts.html GUTS> special cases based on the
% full model (TK and damage dynamics separate): SD, IT and mixed (or GUTS
% proper). 
% 
% *This script:* data set for cadmium in _Daphnia magna_ at 10C from:
% Heugens et al. (2003), <http://dx.doi.org/10.1021/es0264347>. Also
% analysed in: Jager et al (2006),
% <http://dx.doi.org/10.1007/s10646-006-0060-x>. This set has survival data
% and body residues (at a different exposure concentration), which are
% fitted simultaneously. This illustrates how to fit a model without
% explicit damage module when body residues are available.
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
% concentrations) in columns. 

% observed number of survivors, time in days, conc in mg/L
DATA{1} = [ -1 0.00010  0.56  0.82    1.1   1.7   2.2
            0       15    15    15     15    15    15
      0.22917       15    15    15     15    15    15
      0.86458       15    15    15     15    15    15
       1.2292       15    15    15     15    15    15
       1.8576       15    15    15     15    15     9
       2.1736       15    15    15     15    10     6
       2.8854       15    15    15     11     6     3
       3.1979       15    15    13      9     3     3
       3.8646       15    15    13      8     2     0
       4.0417       15    15    13      8     2     0];
 
% actual internal concentrations in mg/kg, time in days
DATA{2} = [ 1    0.099 0.099 0.099
            0    2.907 5 NaN
     0.083333    3.125 6.9277 7.3718
      0.20833    5.5195 14.706 6.7901
      0.33333    8 14 10
            1    16.667 24.324 18.085
        1.875    27.632 25 31.25];   

% scaled damage; there are no data for this stage, but dummy data are added
% to make sure that the damage prediction is plotted by calc_and_plot
% (without having to plot ALL model curves in ALL sub-plots)
DATA{3} = [1 DATA{1}(1,2:end)
           0  1e-6     1e-6    1e-6    1e-6   1e-6   1e-6]; 
       
glo.wts = [1 1 0];  % set zero weight to fake data for scaled damage

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios.

X0mat(1,:) = [DATA{1}(1,2:end) DATA{2}(1,2)]; % scenarios (concentrations or identifiers)
X0mat(2,:) = 1; % initial survival probability
X0mat(3,:) = 0; % initial internal concentration (overwritten by parameter Ci0!
X0mat(4,:) = 0; % initial scaled damage

%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. 

% global parameters for GUTS purposes
glo.sel  = 1; % select death mechanism: 1) SD 2) IT 3) mixed
glo.locS = 1; % location of survival probability in the state variable list
glo.locC = 2; % location of internal concentration in the state variable list
glo.locD = 3; % location of scaled damage in the state variable list
glo.fastrep = 1; % set to 1 to assume fast damage repair (death is driven by Ci)

% start values for SD
% syntax: par.name = [startvalue fit(0/1) minval maxval optional:log/normal scale (0/1)];
par.ke    = [0.7  1 1e-3  10 0]; % elimination rate constant (d-1)
par.kr    = [1.3  1 1e-6  10 0]; % damage repair rate constant (d-1)
par.Kiw   = [370  1    0 1e6 1]; % bioconcentration factor (L/kg)
par.mi    = [200  1    0 1e6 1]; % median threshold for survival (mg/kg)
par.hb    = [0.01 1 1e-6 1e6 1]; % background hazard rate (1/d)
par.bi    = [4e-3 1 1e-6 1e3 0]; % killing rate (kg/mg/d) (SD and mixed)
par.Fs    = [2    1    1 100 1]; % fraction spread of threshold distribution (-) (IT and mixed)
par.Ci0   = [4.1  1    0 100 1]; % initial body residue (mg/kg), overrides values in X0mat!
 
switch glo.sel % make sure that right parameters are fitted
    case 1 % SD
        par.Fs(2) = 0; % never fit the NEC spread!
    case 2 % IT
        par.bi(2) = 0; % never fit the killing rate!
    case 3 % mixed
        % do nothing: fit all parameters
end
if glo.fastrep == 1 % when damage repair is fast ...
    par.kr(2) = 0;  % never fit the repair rate!
end

%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% constructed, based on the data set.

% specify the y-axis labels for each state variable
glo.ylab{1} = 'survival probability';
glo.ylab{2} = 'internal concentration (mg/kg)';
glo.ylab{3} = 'scaled damage (mg/kg)';
% specify the x-axis label (same for all states)
glo.xlab    = 'time (days)';
glo.leglab1 = 'conc. '; % legend label before the 'scenario' number
glo.leglab2 = 'mg/L'; % legend label after the 'scenario' number

prelim_checks % script to perform some preliminary checks and set things up
% Note: prelim_checks also fills all the options (opt_...) with defauls, so
% modify options after this call, if needed.

%% Calculations and plotting
% Here, the function is called that will do the calculation and the
% plotting. The files in this directory always apply the analytical
% solution for damage in simplefun.

opt_optim.fit  = 1; % fit the parameters (1), or don't (0)
opt_optim.it   = 1; % show iterations of the optimisation (1, default) or not (0)
opt_plot.sho   = 0; % set to 1 to show all scenarios, 0 to only plot model for scenarios with data
opt_plot.annot = 2; % extra subplot in multiplot for fits: 1) box with parameter estimates, 2) overall legend

% optimise and plot (fitted parameters in par_out)
par_out = calc_optim(par,opt_optim); % start the optimisation
calc_and_plot(par_out,opt_plot); % calculate model lines and plot them

opt_tktd.obspred = 0; % plot predicted-observed plots (1) or not (0)
opt_tktd.max_exp = 0; % set to 1 to maximise exposure/damage plots on exposure rather than damage
plot_tktd(par_out,opt_tktd,[]); % leaving opt_conf empty suppresses all CIs for these plots

% Below, several options for post analyses are provided.

return % stop here, you can run the next sections afterwards if you like

%% Profiling the likelihood
% By profiling you make robust confidence intervals for one or more of your
% parameters. 
 
opt_prof.detail   = 2;  % detailed (1) or a coarse (2) calculation
opt_prof.subopt   = 10; % number of sub-optimisations to perform to increase robustness
opt_prof.brkprof  = 2; % when a better optimum is located, stop (1) or automatically refit (2)

% % UNCOMMENT FOLLOWING LINE(S) TO CALCULATE
% calc_proflik(par_out,'all',opt_prof,opt_optim);  % calculate a profile
% % Enter single parameter names, a cell array of names, or 'all' to profile
% % all fitted parameters.

%% Likelihood region
% Make intervals on model predictions by using a sample of parameter sets
% taken from the joint likelihood-based conf. region. Sub-optimisations set
% to zero for this demo; in general, it is a good idea to take 10
% sub-optimisations for robustness. Furthermore, the number of elements in
% the sample is now very small, compared to the number of fitted
% parameters. However, doing it properly does require a lot of calculation
% time.

opt_prof.detail  = 1;  % detailed (1) or a coarse (2) calculation
opt_prof.subopt  = 0; % number of sub-optimisations to perform to increase robustness
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
% Options for confidence bounds on model curves can be set using opt_conf
% (see prelim_checks). The plot_tktd function makes multiplots for the
% survival data, which are more readable when plotting with various
% intervals.

opt_conf.type    = 2; % make intervals from 1) slice sampler, 2)likelihood region, 3) parspace explorer
opt_conf.lim_set = 2; % for lik-region sample: use limited set of n_lim points (1) or outer hull (2) to create CIs
opt_tktd.max_exp = 0; % set to 1 to maximise exposure/damage plots on exposure rather than damage
opt_tktd.obspred = 1; % plot predicted-observed plots (1) or not (0)

plot_tktd(par_out,opt_tktd,opt_conf); 

%% Calculate LCx versus time
% Here, the LCx (by default the LC50) is calculated at several time points.
% LCx values are also printed on screen. If a sample from parameter space
% is available (e.g., from the slice sampler or the likelihood region), it
% can be used to calculate confidence bounds. 
% 
% Options for LCx (with confidence bounds) can be set using opt_lcx_lim
% (see prelim_checks). Note that opt_conf.type=-1 skips CIs.

opt_conf.type    = 2; % make intervals from 1) slice sampler, 2)likelihood region, 3) parspace explorer
opt_conf.lim_set = 2; % for lik-region sample: use limited set of n_lim points (1) or outer hull (2) to create CIs

% This is the general method as the fast methods won't work for the full
% model (or at least: won't be faster than the general method).
opt_ecx.Feff      = [0.50]; % effect levels (>0 en <1), x/100 in ECx
opt_ecx.notitle   = 1; % set to 1 to suppress titles above ECx plots
Tend = [1:8]; % times at which to calculate LCx, relative to control

% In this case, the initial concentration Ci0 is also used in the LCx calculation. 
calc_ecx(par_out,Tend,opt_ecx,opt_conf); % general method for ECx values

% If you do not want that, you can set this parameter to zero in the LCx calculation.
opt_ecx.setzero  = {'Ci0'}; % parameter names (as string array) for extra paramaters to be set to zero
calc_ecx(par_out,Tend,opt_ecx,opt_conf); % general method for ECx values

%% Calculate LPx
% We can use the function calc_lpx_lim to calculate an exposure
% multiplication factor: with which factor do we need to multiply an
% exposure scenario to obtain x% effect at the end of the scenario? Here,
% done for 10% effect at the end of the exposure scenario. Note that
% opt_conf.type=-1 skips CIs.
% 
% We have to be a bit careful in this case, since the initial concentration
% is not zero. The fast method of calc_lpx_lim should probably not be used
% as it scales damage with LPx. The general method can be used, either
% including or excluding the initial concentration. The choice will depend
% on the nature of that initial concentration (e.g., whether it partakes in
% TK and can be eliminated or not, and whether it will also occur for
% animals in the field).

opt_conf.type    = 2; % make intervals from 1) slice sampler, 2)likelihood region, 3) parspace explorer
opt_conf.lim_set = 2;  % for lik-region sample: use limited set of n_lim points (1) or outer hull (2) to create CIs
opt_ecx.Feff     = [0.10]; % effect levels (>0 en <1), x/100 in ECx
opt_ecx.notitle  = 1; % set to 1 to suppress titles above ECx plots
opt_ecx.setzero  = {'Ci0'}; % parameter names (as string array) for extra paramaters to be set to zero
Tend  = []; % time at which to calculate LPx (empty is end of profile)

calc_epx(par_out,'profile_monit.txt',Tend,opt_ecx,opt_conf,opt_tktd); % general method for EPx values

