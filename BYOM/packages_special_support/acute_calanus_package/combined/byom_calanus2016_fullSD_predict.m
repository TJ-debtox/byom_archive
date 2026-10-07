%% BYOM, byom_calanus2016_fullSD_predict.m GUTS
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
% *The model:* fitting survival data and body residues simultaneously with
% the <http://www.debtox.info/about_guts.html GUTS> special cases based on
% the internal concentration: SD, IT and mixed. The internal concentration
% in structure is directly linked to the effect, without damage stage.
%
% *This script:* C2-naphthalene in _Calanus finmarchicus_. Fitting both the
% survival data and the two-compartment TK model on the body residues.
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

DATA{1} = [0];
DATA{2} = [0];
DATA{3} = [0];
DATA{4} = [0];

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios.

X0mat(1,:) = [1.2]; % scenarios (concentrations)
X0mat(2,:) = 1; % initial survival probability
X0mat(3,:) = 0; % initial internal concentration in structure
X0mat(4,:) = 0; % initial internal concentration in lipids
X0mat(5,:) = 0; % initial internal concentration in total animal

%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. 

% global parameter
glo.FLS  = 0.2; % VL/VS
glo.Tdep = 10; % depuration time for simulations
glo.sel  = 1; % select death mechanism: 1) SD 2) IT
% globals used for GUTS purposes
glo.locD = 2; % location of the dose metric in the state variable list
glo.locS = 1; % location of survival probability in the state variable list

% syntax: par.name = [startvalue fit(0/1) minval maxval optional:log/normal scale (0/1)];
par.keS    = [   0.55301  1       0.01         10 1]; % elimination rate constant for structure (d-1)
par.kuS    = [    1442.2  1       0.01      1e+06 1]; % uptake rate constant for structure (L/L/d)
par.keL    = [      0.01  1       0.01         10 1]; % elimination rate constant for lipids (d-1)
par.kuL    = [    0.8248  1       0.01      1e+06 1]; % uptake rate constant for lipids (L/L/d)
par.mS     = [    7942.6  1          0      1e+06 1]; % median threshold for survival (umol/L)
par.hb     = [ 0.0045103  1          0      1e+06 1]; % background hazard rate (1/d)
par.bS     = [0.00025968  1          0       1000 1]; % killing rate (L/umol/d)
par.Fs     = [       1.5  0          1        100 1]; % fraction spread of NEC distribution (-) (IT and mixed)

switch glo.sel % make sure that right parameters are fitted
    case 1 % SD
        par.Fs(2) = 0; % never fit the NEC spread!
    case 2 % IT
        par.bS(2) = 0; % never fit the killing rate!
end

%% Time vector and labels for plots
% Specify what to calculate and what to plot.

glo.t   = linspace(0,20,100); % time vector for the model curves in days

% specify the y-axis labels for each state variable
glo.ylab{1} = 'survival probability';
glo.ylab{2} = ['internal conc., structure (',char(181),'mol/L)'];
glo.ylab{3} = ['internal conc., lipid sac (',char(181),'mol/L)'];
glo.ylab{4} = ['internal conc., total (',char(181),'mol/kg)'];

% specify the x-axis label (same for all states)
glo.xlab    = 'time (days)';
glo.leglab1 = ''; % legend label before the 'scenario' number
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

glo.stiff    = [0 1]; % ODE solver 0) ode45 (standard), 1) ode113 (moderately stiff), 2) ode15s (stiff)
% Seems that standard solver performs good enough, with tight tolerances.
% Second argument is for normally tight (1), tighter (2), or very tight (3)
% tolerances.
glo.break_time = 0; % break time vector up for ODE solver (1) or don't (0)

opt_optim.fit  = 0; % fit the parameters (1), or don't (0)
opt_optim.it   = 0; % show iterations of the simplex optimisation (1, default) or not (0)
opt_plot.sho   = 1; % set to 1 to show all scenarios, 0 to only plot model for scenarios with data
opt_plot.bw    = 1; % plot in black and white
opt_plot.cn    = 0; % if set to 1, connect model line to points (only for bw=1)

% plot (fitted parameters in par_out)
par_out = par; % no optimisation, just use the start values for the simulations

opt_plot.statsup = [1 2 3]; % vector with states to suppress in plotting fits
calc_and_plot(par_out,opt_plot); % calculate model lines and plot them

ylim([0 6000])
str = ['$V_L/V_S = $',num2str(glo.FLS)];
text(gca,10,3000,str,'Interpreter','latex','FontSize',18)

opt_plot.statsup = [1 3 4]; % vector with states to suppress in plotting fits
calc_and_plot(par_out,opt_plot); % calculate model lines and plot them

ylim([0 3500])
str = ['$V_L/V_S = $',num2str(glo.FLS)];
text(gca,10,3000,str,'Interpreter','latex','FontSize',18)

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
Tend = [1:0.2:3 3.5 4 5 6 8 10]; % times at which to calculate LCx, relative to control

calc_ecx(par_out,Tend,opt_ecx,opt_conf); % general method for ECx values

ylim([0 15])
str = ['$V_L/V_S = $',num2str(glo.FLS)];
text(gca,4,12,str,'Interpreter','latex','FontSize',18)