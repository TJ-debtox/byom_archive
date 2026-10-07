%% BYOM, byom_guts_LEV.m, full GUTS model
%
% * Author: Tjalling Jager 
% * Date: November 2023
% * Web support: <http://www.debtox.info/byom.html>
%
% BYOM is a General framework for simulating model systems. The files in
% this directory use an analytical solution only, and therefore
% <derivatives.html derivatives.m> will be missing.
%
% *The model:* fitting survival data with the
% <http://www.debtox.info/about_guts.html GUTS> special cases based on the
% full model (TK and damage dynamics separate): SD, IT and mixed (or GUTS
% proper). The standard full GUTS model package was adapted slightly to
% allow the double-scaled model version, and to easily switch from the
% script file between the full model and the reduced model (using
% glo.fastkin). 
% 
% *This script:* data set for LEV in _Daphnia magna_ from: Tolosi and De
% Liguoro (2021), <https://doi.org/10.1016/j.ecoenv.2021.112778>. This
% script allows to redo the calibrations from my manuscript on delayed
% toxicity (in prep.). The data points were extracted from the original
% publication using <https://github.com/jornbr/plotreader/wiki PlotReader>,
% so there may be small errors (particularly where treatments overlap in
% low-level response).
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
glo.saveplt = 3; % save all plots as (1) Matlab figures, (2) JPEG file or (3) PDF (see all_options.txt)

%% The data set
% Data are entered in matrix form, time in rows, scenarios (exposure
% concentrations) in columns. 

% observed number of survivors, time in days, conc in mg/L
DATA{1} = [-1	0	0.7	1.2	2.1	3.8	6.9	12.3	22.2	40
0	20	20	20	20	20	20	20	20	20
1	20	20	20	20	20	20	20	20	20
2	20	20	20	20	20	20	20	19	20
3	20	20	20	20	20	20	20	19	20
5	20	20	20	20	20	19	18	16	11
7	20	20	20	20	20	19	17	11	7
9	20	20	20	20	20	19	16	9	4
12	20	20	20	20	19	18	14	5	1];
   
% DATA{1} = DATA{1}(1:4,:); % uncomment to restrict data set to first 2 days!

% Define the exposure scenarios as a series of pulses with a new constant
% exposure in each time period. This allow use of an analytical solution
% for the TK/damage dynamics in simplefun.
Cw = [0 0.7	1.2	2.1	3.8	6.9	12.3	22.2	40
    0 0.7	1.2	2.1	3.8	6.9	12.3	22.2	40
    2 0       0   0   0   0    0       0     0];
 
make_scen(2,Cw); % create the globals to define the forcing function

DATA{2} = []; % internal concentrations are not provided

% scaled damage; there are no data for this stage, but dummy data are added
% to make sure that the damage prediction can be plotted by calc_and_plot
% (without having to plot ALL model curves in ALL sub-plots)
DATA{3} = [1 DATA{1}(1,2:end)
           0  1e-6     1e-6    1e-6    1e-6   1e-6   1e-6    1e-6   1e-6   1e-6]; 
       
glo.wts = [1 1 0];  % set zero weight to fake data for scaled damage

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios.

X0mat(1,:) = [DATA{1}(1,2:end)]; % scenarios (concentrations or identifiers)
X0mat(2,:) = 1; % initial survival probability
X0mat(3,:) = 0; % initial (scaled) internal concentration
X0mat(4,:) = 0; % initial scaled damage

%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. 

% global parameters for GUTS purposes
glo.sel  = 2; % select death mechanism: 1) SD 2) IT 3) mixed
glo.locS = 1; % location of survival probability in the state variable list
glo.locC = 2; % location of internal concentration in the state variable list
glo.locD = 3; % location of scaled damage in the state variable list
glo.fastkin = 0; % set to 1 to assume fast toxicokinetics (death is still driven by Dw)

% Create a base-name for the MAT files and figures that are produced by
% this run. The name is based on the settings, so we can keep the apart.
fname_str = [mfilename,'_',num2str(DATA{1}(end,1)),'d'];
switch glo.sel
    case 1
        fname_str = [fname_str,'_SD'];
    case 2
        fname_str = [fname_str,'_IT'];
    case 3
        fname_str = [fname_str,'_MIX'];
end
if glo.fastkin == 1
    fname_str = [fname_str,'_1cmp'];
else
    fname_str = [fname_str,'_2cmp'];
end
glo.basenm  = fname_str; % use this as filename for for the plots and MAT file

save([glo.basenm,'_DATA'],'DATA','W') % save MAT file with data set 
% Saving the data set is handy to allow for a simple reconstruction of the
% calibrations, without needing to define the data again, in the same way.
% Note that if you use DATAx and Wx, you need to save them as well!

if glo.sel == 1 % start values and bounds for SD
    % syntax: par.name = [startvalue fit(0/1) minval maxval optional:log/normal scale (0/1)];
    par.kd    = [0.1     1 1e-3  100 0]; % dominant rate constant (d-1)
    par.ke    = [0.35    1 1e-3   10 0]; % elimination rate constant (d-1)
    par.kr    = [0.001   1 1e-3    1 0]; % damage repair rate constant (d-1)
    % NOTE: we can reverse ke and kr and get a good fit, so it unclear which of
    % the two processes is really slow. Here, I force ke>kr in call_deri.
    par.mw    = [0.0030  1 1e-6   10 0]; % median threshold for survival (mg/L)
    par.hb    = [1e-3    1 1e-4 1e-1 0]; % background hazard rate (1/d)
    par.bw    = [8.0     1 1e-2  1e4 0]; % killing rate (L/mg/d) (SD and mixed)
    par.Fs    = [2       1    1   20 1]; % fraction spread of threshold distribution (-) (IT and mixed)

else % start values and bounds for IT
    % syntax: par.name = [startvalue fit(0/1) minval maxval optional:log/normal scale (0/1)];
    par.kd    = [0.1     1 1e-3  100 0]; % dominant rate constant (d-1)
    par.ke    = [0.35    1 1e-3    1 0]; % elimination rate constant (d-1)
    par.kr    = [0.001   1 1e-3 1e-1 0]; % damage repair rate constant (d-1)
    % NOTE: we can reverse ke and kr and get a good fit, so it unclear which of
    % the two processes is really slow. Here, I force ke>kr in call_deri.
    if glo.fastkin == 1
        par.mw    = [0.0030  1 1e-5  100 0]; % median threshold for survival (mg/L)
        par.hb    = [1e-3    1 1e-4 1e-1 0]; % background hazard rate (1/d)
        par.bw    = [8.0     1 1e-2  1e4 0]; % killing rate (L/mg/d) (SD and mixed)
        par.Fs    = [2       1    1   30 1]; % fraction spread of threshold distribution (-) (IT and mixed)
    else % for the full model, other bounds are more efficient
        par.mw    = [0.0030  1 1e-5  5 0]; % median threshold for survival (mg/L)
        par.hb    = [1e-3    1 1e-4 1e-1 0]; % background hazard rate (1/d)
        par.bw    = [8.0     1 1e-2  1e4 0]; % killing rate (L/mg/d) (SD and mixed)
        par.Fs    = [2       1    1   10 1]; % fraction spread of threshold distribution (-) (IT and mixed)
    end

end

switch glo.sel % make sure that right parameters are fitted
    case 1 % SD
        par.Fs([1 2]) = [1 0]; % never fit the threshold spread and set to 1
    case 2 % IT
        par.bw([1 2]) = [par.bw(4) 0]; % never fit the killing rate
    case 3 % mixed
        % do nothing: fit all parameters
end
if glo.fastkin == 1 % when toxicokinetics is fast ...
    par.kr([1 2]) = [par.kr(4) 0];  % never fit the repair rate, and set to max of range
    par.ke([1 2]) = [par.ke(4) 0];  % never fit the elimination rate, and set to max of range
else 
    par.kd([1 2]) = [par.kd(4) 0];  % never fit the dominant rate, and set to max of range
end

%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% constructed, based on the data set.

% specify the y-axis labels for each state variable
glo.ylab{1} = 'survival prob. (-)';
glo.ylab{2} = 'scaled int. conc. (mg/kg)';
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

opt_optim.type   = 4; % optimisation method 1) simplex, 4 parameter-space explorer
opt_optim.fit    = 1; % fit the parameters (1), or don't (0)
opt_plot.bw      = 1; % if set to 1, plots in black and white with different plot symbols
opt_plot.annot   = 2; % annotations in multiplot for fits: 1) box with parameter estimates 2) single legend
opt_plot.statsup = [glo.locC glo.locD]; % vector with states to suppress in plotting fits

opt_optim.ps_plots = 0; % when set to 1, makes intermediate plots to monitor progress of parameter-space explorer
opt_optim.ps_profs = 1; % when set to 1, makes profiles and additional sampling for parameter-space explorer
opt_optim.ps_rough = 0; % set to 1 for rough settings of parameter-space explorer, 0 for settings as in openGUTS
opt_optim.ps_saved = 0; % use saved set for parameter-space explorer (1) or not (0);
opt_optim.ps_dupl  = 0; % set to 1 to remove duplicates from sample in calc_parspace (slow for rough=0!)
opt_optim.ps_notitle = 1; % set to 1 to suppress plotting title on parameter-space plot

% optimise and plot (fitted parameters in par_out)
par_out = calc_optim(par,opt_optim); % start the optimisation
Xall = calc_and_plot(par_out,opt_plot); % calculate model lines and plot them

% print_par(par_out)

opt_conf.type    = 3; % make intervals from 1) slice sampler, 2)likelihood region, 3) parspace explorer
opt_conf.lim_set = 2; % for lik-region sample: use limited set of n_lim points (1) or outer hull (2) to create CIs
opt_tktd.max_exp = 0; % set to 1 to maximise exposure/damage plots on exposure rather than damage
opt_tktd.obspred = 0; % plot predicted-observed plots (1) or not (0)
opt_tktd.flip    = 0; % set to 1 to flip row 1 and 2 of the plot around
opt_tktd.statsup = [glo.locC glo.locD]; % states to suppress from the plots (e.g., locS)
% opt_tktd.statsup = []; % states to suppress from the plots (e.g., locS)
opt_tktd.notitle  = 1; % set to 1 to suppress titles above plots

h_coll = plot_tktd(par_out,opt_tktd,opt_conf); % leaving opt_conf empty suppresses all CIs for these plots

% resize the multipanel plot for survival and all treatments to a more
% practical size for in the supporting information
figh = h_coll{1};
figh.Position(3) = 2000;
figh.Position(4) = 250;

% overwrite the previous plot made by plot_tktd with the resized one
savenm = ['tktd_plot_',glo.basenm];%
save_plot(figh,savenm);

