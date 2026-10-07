%% BYOM, byom_guts_ttd_test.m
%
% * Author: Tjalling Jager
% * Date: September 2023
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
% *This script:* Demonstrating the use of data in a time-to-death format
% rather than the more typical survivor counts over time.
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
% * -4 for time-to-death data
% * -1 for multinomial likelihood (for survival data)
% * 0  for log-transform the data, then normal likelihood
% * 0.5 for square-root transform the data, then normal likelihood
% * 1  for no transformation of the data, then normal likelihood

% Time-to-death data are entered in a matrix, where each cell represents
% the time of death for an individual. A negative value is used for when
% the animal did not die: it should be set to the last time point at which
% we know it was still alive. In the hypothetical example below, the test
% duration is 4 days, so -4 implies that an individual was alive at the end
% of the test. A -2.2 implies that an individual disappeared or was taken
% out after 2.2 days. If the various treatments do not have equal numbers
% of individuals, NaN's can be used to fill up the matrix. First row are
% the treatments (as is BYOM standard). First column (usually time) is now
% arbitary: it is not used, but it is probably a good idea to make the
% values unique and increasing.
% 
% NOTE: when plotting the data, plot_tktd will translate the TTD data into
% survivors over time. The number of time points is set in the option
% opt_tktd.ttd_stp.

% Hypothetical data set to demonstrate how to use TTD data sets
TTD = [-4 0 18 100
    1   -4	3.1	1.2
    2   -4	2.5	3.4
    3   3.2	3.9	0.7
    4   -4	-4 -2.2
    5   -4	NaN	1.5
    6   3.6	2.9	-1.1
    7   -4	1.4	2.8];

DATA{1} = TTD;

% % Alternatively, we can translate the TTD matrix into a regular matrix (and
% % a weights matrix for the missing/removed individuals). This gives a very
% % similar fit, but a different absolute value for the likelihood function.
% [S,w]  = conv_ttd_matrix(TTD,100); % this forces the continuous observations into 19 regular time intervals
% DATA{1} = S;
% W{1}    = w;

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
par.kd = [0.001  1 1e-3  100 0];   % dominant rate constant, d-1
par.mw = [0.006  1   0  1e6 1];   % median threshold for survival (ug/L)
par.hb = [0.07   1    0    1 1];   % background hazard rate (1/d)
par.bw = [4.5    1 1e-6  1e6 0];   % killing rate (L/ug/d) (SD only)
par.Fs = [3      1    1  100 1];   % fraction spread of threshold distribution (IT only)

% % start values for IT
% % syntax: par.name = [startvalue fit(0/1) minval maxval optional:log/normal scale (0/1)];
% par.kd    = [0.001 1 1e-3 100 0]; % dominant rate constant (d-1)
% par.mw    = [0.2   1    0 1e6 1]; % median threshold for survival (ug/L)
% par.hb    = [0.1   1    0   1 1]; % background hazard rate (1/d)
% par.bw    = [0.04  1 1e-6 1e6 0]; % killing rate (L/ug/d) (SD and mixed)
% par.Fs    = [5     1    1 100 1]; % fraction spread of NEC distribution (-) (IT and mixed)

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
% prelim_checks.m). Options for the optimsation routine can be set using
% opt_optim. The files in this directory always apply the analytical
% solution for damage in simplefun.
%
% For the demo, the iterations were turned off (opt_optim.it = 0). In this
% case, the default position of the legend is obscuring part of the data
% and model lines. In <byom_guts_extra.html byom_guts_extra.m> an option is
% demonstrated to place the legend in a separate subplot.

% par = start_vals_guts(par); % experimental start-value finder; use at your own risk
% % Note: start_vals will now overwrite the parameter structure par!

opt_optim.fit = 1; % fit the parameters (1), or don't (0)
opt_optim.it  = 1; % show iterations of the optimisation (1, default) or not (0)

% optimise and plot (fitted parameters in par_out)
par_out = calc_optim(par,opt_optim); % start the optimisation

opt_plot.annot   = 1; % extra subplot in multiplot for fits: 1) box with parameter estimates, 2) overall legend
% calc_and_plot(par_out,opt_plot); % calculate model lines and plot them
opt_tktd.obspred = 0; % plot predicted-observed plots (1) or not (0)
plot_tktd(par_out,opt_tktd,[]); % make more detailed standard TKTD plots (without CIs)
