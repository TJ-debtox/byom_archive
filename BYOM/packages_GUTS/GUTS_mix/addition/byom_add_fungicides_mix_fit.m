%% BYOM, byom_add_fungicides_mix_fit.m
%
% * Author: Tjalling Jager
% * Date: December 2021
% * Web support: <http://www.debtox.info/byom.html>
%
% BYOM is a General framework for simulating model systems in terms of
% ordinary differential equations (ODEs). The model itself needs to be
% specified in <derivatives.html derivatives.m> or <simplefun.html
% simplefun.m>, and <call_deri.html call_deri.m> may need to be modified to
% the particular problem as well. The files in the engine directory are
% needed for fitting and plotting. Results are shown on screen but also
% saved to a log file (|results.out|).
%
% *The model:* fitting binary mixtures for survival data with the
% <http://www.debtox.info/about_guts.html GUTS> special cases based on the
% reduced model (TK and damage dynamics lumped): SD or IT. All
% calculations performed in <simplefun.html simplefun.m>. This is for
% damage addition.
%
% *This script:* Data for prochloraz and triflumizole in Enchytraeus
% crypticus from Bart et al (2021). The files in this directory use the
% analytical solution only! This script is set up for fitting on the entire
% mixture data set.
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

% Treatments are referred to by a two-digit code as in xy, which would mean
% concentration number x for chemical A and concentration number y for
% chemical B. Here, the first chemical is IMI (A) and second THI (B). The
% global variable tab_mix is a translation table for simplefun. In the
% future, I would like to implement this deeper into BYOM, probably into
% make_scen. Note that you can use any coding that you like, as long as
% each code refers to the correct concentration combination (which is laid
% down in tab_mix).
%
% The code below makes use of the option in BYOM to have multiple data sets
% per state, so the single exposures of IMI are a separate data set from
% the single exposures of THI. This means they are plotted in separate
% graphs, but has no consequences for the calculations. For mixture
% survival, you can add DATA{3,1}, or more data sets, with their own unique
% treatment coding.
%  
% The exposure definition is done with make_scen. All possible single
% exposures are defined; the mixtures are all a combination of single
% exposures. In simplefun, the scenario identifier is separated into a
% scenario for chemical A and a scenario for chemical B. So, scenario 42
% implies taking exposure to chemical A from scenario 40, and exposure to
% chemical B from scenario 2. By default, we use a factor of 10 to separate
% the scenario into chemical A and B. This can be modified to 100, by
% changing glo.mix_fact. In that case, IDs 1-99 are for chemical A, and
% 100-9900 (in steps of 100) for chemical B.

glo.mix_fact  = 10; % factor for scenario IDs of the mixture analysis
glo.scen_plot = 0; % don't make plots for the scenarios; it's constant exposure anyway

% Chemical A is prochloraz (mg/L)
CwA = [1 10	 20	 30	 40
       0 10  25  40  60
       2 10  25  40  60];

make_scen(2,CwA);

% Chemical B is triflumizole (mg/L)
CwB = [1  1	  2	  3	  4
       0 10  25  40  60
       2 10  25  40  60];

make_scen(2,CwB);
   
% Make scenarios for the mixture design (they don't overlap with the
% singles here!)
CwAm = [1 50   60 70 80
       0   5 12.5 20 30
       2   5 12.5 20 30];

make_scen(2,CwAm);

CwBm = [1  5    6  7  8
       0   5 12.5 20 30
       2   5 12.5 20 30];

make_scen(2,CwBm);

% Survivor data for prochloraz (A) only        
DATA{1,1} = [-1	  0  10  20  30	 40
              0  24  24  24  24  24
          0.125  24  24  24  24  24
           0.25  24  24  24  24  24
              1  24  24  21  13   0
              2  24  24  17   1   0
              3  24  24   5   0   0
              4  24  24   2   0   0];

% Survivor data for triflumizole (B) only
DATA{2,1} = [-1	  0   1   2   3	  4
              0  24  24  24  24  24
          0.125  24  24  24  23  16
           0.25  24  24  24  17   3
              1  24  24  24  12   0
              2  24  24  24   2   0
              3  24  24  23   2   0
              4  24  24  22   2   0];

% Survivor data for the mixture (here only equiconcentration mixtures). The
% mixture data may also be split up over multiple data sets.
DATA{3,1} = [-1 0 55 66 77 88
              0	24	24	24	24	24
          0.125	24	24	24	23	17
           0.25	24	24	23	23	10
              1	24	24	21	14	1
              2	24	24	12	0	0
              3	24	24	11	0	0
              4	24	22	11	0	0];


%% Create a table with nicer labels for the legends
% Creating a Matlab table in glo.LabelTable will replace the automatically
% generated labels (from glo.leglab1/2 and the scenario identifiers) with a
% dedicated text for each scenario.

% Start with prochloraz single exposure
Scenario = [0 10 20 30 40]';
Label    = {'control';'PRO 10 mg/L';'PRO 25 mg/L';'PRO 40 mg/L';'PRO 60 mg/L'};
% And add the triflumizole single exposures (don't add control again)
Scenario = [Scenario;[1 2 3 4]'];
Label    = [Label;{'TRI 10 mg/L';'TRI 25 mg/L';'TRI 40 mg/L';'TRI 60 mg/L'}];
% And add the mixture exposure (don't add control again)
Scenario = [Scenario;[55 66 77 88]'];
Label    = [Label;{'MIX 5+5 mg/L';'MIX 12.5+12.5 mg/L';'MIX 20+20 mg/L';'MIX 30+30 mg/L'}];

glo.LabelTable = table(Scenario,Label); % create a Matlab table for the labels

%% Call the Matlab GUI open-file element to load MAT file

[~,~,par] = select_pred([1 1 0]); 

% Note that this function, with these settings, loads a pre-saved MAT file.
% The parameter structure <par> is loaded from the selected MAT file. This
% contains starting values for the parameter-space explorer, as generated
% by <combine_mat_files>. It also defines glo.sel.

glo = rmfield(glo,'mat_nm'); % but remove the name of the MAT file as we'll make a new one

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios.

X0mat = glo.LabelTable.Scenario'; % scenarios to be run (treatment codes)
X0mat(2,:) = 1;            % initial survival probability
X0mat(3,:) = 0;            % initial scaled damage
X0mat(4,:) = 0;            % initial scaled damage

%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. 

% global parameters for GUTS purposes
glo.locS = 1; % location of survival in the state-variable vector
glo.locD = [2 3]; % location of damage states in the state-variable vector
% Note that glo.sel is read from the file name by select_pred.

% do fit hb
par.hb = [0.01 1 1e-4    1 1];   % background hazard rate (1/d)
% Note: in this case, it is important to NOT let hb go to zero. There is on
% animal that dies in the first mixture treatment, and the data set can
% only be properly fitted under the assumption that this background
% mortality and NOT related to the chemical(s).

%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% constructed, based on the data set.

% specify the y-axis labels for each state variable
glo.ylab{1} = 'survival probability';
glo.ylab{2} = 'scaled damage PRO (mg/L)';
glo.ylab{3} = 'scaled damage TRI (mg/L)';
% specify the x-axis label (same for all states)
glo.xlab    = 'time (days)';
glo.leglab1 = 'scen. '; % legend label before the 'scenario' number
glo.leglab2 = ''; % legend label after the 'scenario' number

prelim_checks % script to perform some preliminary checks and set things up
% Note: prelim_checks also fills all the options (opt_...) with defauls, so
% modify options after this call, if needed.

% % Below a little trick to NOT fit ANY tox parameter anymore, apart from the
% % interaction factor. That way, you can just fit that one and see if it is
% % significantly different from 0. This must be placed below prelim_checks
% % since that function defines glo2.
% for i = 1:length(glo2.names)
%     if ~strcmp(glo2.names{i},'hb') % hb should stay as is
%         par.(glo2.names{i})(2) = 0;
%     end
% end
% % but DO fit the interaction factor and hb
% par.IAB    = [0  1 -10 10 1]; % interaction factor

%% Calculations and plotting
% Here, the function is called that will do the calculation and the plotting.
% Options for the plotting can be set using opt_plot (see prelim_checks.m).
% Options for the optimisation routine can be set using opt_optim. 

opt_optim.fit    = 1; % fit the parameters (1), or don't (0)
opt_optim.it     = 0; % show iterations of the optimisation (1, default) or not (0)
opt_plot.bw      = 1; % if set to 1, plots in black and white with different plot symbols
opt_plot.annot   = 2; % annotations in multiplot for fits: 1) box with parameter estimates 2) single legend

% The switch fit_tox is used to select whether to fit the control treatment
% (0), the toxicity treatments (1; the control is shown as well), or both
% together (2). A good strategy is to fit the background hazard to the
% control data first (fit_tox=0), copy the best values into the parameter
% matrix above, and then fit the toxicity parameter to the complete data
% set (fit_tox=1). The code below automatically does that.

skip_sg = 1; % set to 1 to skip startgrid completely (use ranges in <par> structure)
% Note: here we do need to skip startgrid as it only works for standard
% GUTS analyses, and we have starting ranges from the saved MAT file.

% ===== FITTING CONTROLS ==================================================
% Simplex fitting works fine for background hazard
opt_optim.type = 1; % optimisation method: 1) default simplex, 4) parspace explorer
fit_tox = [0]; % Fit background hazard rate in data set
opt_plot.statsup = [2 3]; % vector with states to suppress in plotting fits
par_out = automatic_runs_basic(fit_tox,par,-1,skip_sg,opt_optim,opt_plot);
% script to run the calculations and plot, automatically
par = copy_par(par,par_out,1); % copy fitted parameters into par, and keep fit mark in par

% ===== FITTING TOX DATA ==================================================
% Here, we use the parameter-space explorer for fitting the treatments.
% Note that, with fit_tox=1, the tox parameters in par are replaced (when
% fitted) with estimates based on the data set using startgrid_guts
% (called in automatic_runs). You can set skip_sg=0 to use the ranges in
% par (as defined above) instead.

opt_optim.type     = 4; % optimisation method: 1) default simplex, 4) parspace explorer
opt_optim.ps_plots = 0; % when set to 1, makes intermediate plots of parameter space to monitor progress
opt_optim.ps_profs = 1; % when set to 1, makes profiles and additional sampling for parameter-space explorer
opt_optim.ps_rough = 1; % set to 1 for rough settings of parameter-space explorer, 0 for settings as in openGUTS
% Note: skipping profiles and rough settings will speed up the
% optimisation. This is good for exploration of the data. However, always
% best to do the final analysis with full detail.

% create a new basenm such that the death mechanism is added (this makes it
% easier to use the predict script in this folder)
basenm_rem = glo.basenm; % remember basename as we will modify it!
glo.basenm = [basenm_rem,'_sel',sprintf('%d',glo.sel)]; 

fit_tox = [1]; % Fit tox parameters in data set
opt_plot.statsup = []; % vector with states to suppress in plotting fits
par_out = automatic_runs_basic(fit_tox,par,-1,skip_sg,opt_optim,opt_plot);


% % optimise and plot (fitted parameters in par_out)
% par_out = calc_optim(par,opt_optim); % start the optimisation
% calc_and_plot(par_out,opt_plot); % calculate model lines and plot them

% print_par(par_out) % display fitted parameters on screen in formatted manner to they can be copied into the code

% % If you want to make dedicated plots yourself (e.g., a predicted vs
% % observed plot), note that calc_and_plot can also return the numerical
% % output for the state variables. See example below. With Xall, you can
% % make your own plots.
% glo.t = [0 1 2]; % only calculate the time points of the data set
% Xall = calc_and_plot(par_out,opt_plot); % calculate model lines and plot them

%% Plot results with confidence intervals
% The following code can be used to make plots with confidence intervals.
% Options for confidence bounds on model curves can be set using opt_conf
% (see prelim_checks). The plot_tktd function makes multiplots for the
% survival data, which are more readable when plotting with various
% intervals.

opt_conf.type    = 3; % make intervals from 1) slice sampler, 2)likelihood region, 3) parspace explorer
opt_conf.lim_set = 2; % use limited set of n_lim points (1) or outer hull (2, not for Bayes) to create CIs
opt_tktd.notitle = 1; % set to 1 to suppress titles above plots (becomes messy for GUTS)
opt_tktd.obspred = 2; % plot predicted-observed plots (1) or not (0), (2) makes 1 plot for multiple data sets

% Calculate and plot dedicated TKTD plots. These plots are more readable
% when plotting CIs.
plot_tktd([],opt_tktd,opt_conf); 
% Note that first input is left empty: this is for the parameter structure,
% which is then obtained from the saved sample.
% leave the options opt_conf empty to suppress all CIs for these plots


