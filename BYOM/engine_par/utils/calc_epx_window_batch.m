function calc_epx_window_batch(par_plot,Twin,opt_ecx,opt_conf)

% Usage: calc_epx_window_batch(par_plot,Twin,opt_ecx,opt_conf)
%           PART OF ENGINE_PAR
% Batch calculations for EPx values due to an exposure profile, with moving
% time window. This functions uses the Matlab GUI element to select files,
% and calls calc_epx_window_sub (which is a sub-function based on a
% stripped version of calc_epx_window) repeatedly for the actual
% calculations. If opt_conf is provided, and CIs are requested, CIs are
% calculated on the lowest EPx for each profile, each trait, and each
% effect level x. This is only done for the window in which the lowest EPx
% is found (for that trait and that effect level). This is much faster than
% calculating CIs for each window (which is what calc_epx_window does).
% 
% <par_plot>   parameter structure for the best-fit curve; if left empty the
%              structure from the saved sample is used
% <Twin>       length of time window (days)
% <opt_ecx>    options structure for ECx and EPx calculations
% <opt_conf>   options structure for making confidence intervals (needed if
%              par_plot is to be read from file)
% 
% Currently, this function does not provide any output, outer than to
% screen (and to the log file).
% 
% Author     : Tjalling Jager 
% Date       : July 2023
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2024, Tjalling Jager.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE-MIT.txt file in the root directory of BYOM. 

global glo

mf_crit   = opt_ecx.mf_crit;  % MF trigger for flagging potentially critical profiles
calc_int  = opt_ecx.calc_int; % integrate survival and repro into 1) RGR, or 2) survival-corrected repro (experimental!)
id_sel    = opt_ecx.id_sel;   % scenario to use from X0mat, scenario identifier, flag for ECx to use scenarios rather than concentrations
Feff      = opt_ecx.Feff;     % effect levels (>0 en <1), x/100 in ECx/EPx
par_read  = opt_ecx.par_read; % when set to 1 read parameters from saved set, but do NOT make CIs
rob_win   = opt_ecx.rob_win;  % set to 1 to use robust EPx calculation for moving time windows, rather than with fzero
saveall   = opt_ecx.saveall;   % set to 1 to save all output (EPx for each element of the sample) in separate folder
                               % NOTE: if multiple profiles are selected, this option is a bad idea!
max_mf    = opt_ecx.max_mf;   % maximum allowed multiplication factor

if isempty(opt_conf)
    type_conf = 0;
else
    type_conf = opt_conf.type; % use values from slice sampler (1), likelihood region(2) to make intervals
    type_conf = max(0,type_conf); % if someone uses -1, set it to zero
end

% see which states are there; initialise with silly high values so they
% won't be set unless they exist (and switch-case runs well)
locS  = 100;
locL  = 101;
locR  = 102;
loc_a = 103;
if isfield(glo,'locS') % then we have a state of survival
    locS = glo.locS; % collect the location
end
if isfield(glo,'locL') % then we have a state of body length
    locL = glo.locL; % collect the location
end
if isfield(glo,'locR') % then we have a state of reproduction
    locR = glo.locR; % collect the location
end
if isfield(glo,'loc_a') % then we have a state of active
    loc_a = glo.loc_a; % collect the locations
end
loc_all = [locS locL locR loc_a];

%% Select and load profiles from text
% These may be located in a different location than the current working
% directory. Therefore, filepath is also used.

diary(glo.diary) % collect any warnings in the diary "results.out"
warning('off','backtrace')
warn_init = [0 0];
if ~isempty(glo.names_sep)
    warning('You are using separate parameters per data set (with glo.names_sep)')
    warning(['Parameters for data set ',num2str(1+floor(id_sel(2)/100)),' will be used for effects'])
    disp(' ')
    warn_init(1) = 1;
end
if rob_win == 0 && isfield(glo,'moa') && isfield(glo,'feedb')
    if glo.feedb(2) == 1 && any(glo.moa(1:3)==1)
        warning('Calculating standard EPx windows is NOT recommended for this combination of pMoA and feedbacks!')
        warning('There is a possibility of (limited) lower EPx values than those reported.')
        warning('It is advisable to use the robust settings with opt_ecx.rob_win = 1.')
        disp(' ')
        warn_init(2) = 1;
    end
end
warning('on','backtrace')
diary off

% Use Matlab GUI element to select files for loading
[filename,filepath] = uigetfile('input_data/*.txt','Select file(s) with exposure profiles for batch analysis','MultiSelect','on'); % use Matlab GUI to select file(s)
if ~iscell(filename) && numel(filename) == 1 && filename(1) == 0 % if cancel is pressed ...
    return % simply stop
end
if ~iscell(filename) % the GUI makes it a cell array only when multiple files are selected ...
    filename = {filename}; % but we also want a cell array if it's only one
end

if length(filename) > 1 && saveall == 1 
    % if multiple profiles are selected, this option is a bad idea because
    % all output is saved into the same folder. It is thus unclear which
    % output belongs to which profile, and output may overwrite earlier
    % files. So set to 0!
    warning('off','backtrace')
    warning('Using the saveall option with multiple exposure files is not a good idea, so it is turned off!')
    warning('on','backtrace')
    opt_ecx.saveall = 0; % set to 1 to save all output (EPx for each element of the sample) in separate folder
end

% If we need CIs, load the best parameter set and the random sample from
% file (that way, we can feed rnd to calc_epx, so that function does not
% have to load it again for every window).
if type_conf > 0 || isempty(par_plot) % also if par_plot is not provided
    if isempty(opt_conf)
        error('opt_conf needs to be provided to read random sample and/or parameter structure from MAT file')
    end
    [rnd,par] = load_rnd(opt_conf); % loading and selecting the sample is handled in a separate function
    if numel(rnd) <= 1 % that means that no (useful) sample was found
        type_conf = 0; % no need to produce an error, just do analysis without CI
    end
    if isempty(par_plot) % if no par structure was entered in this function ...
        par_plot = par; % simply use the one from the sample file
        disp('The parameter structure is read from the MAT file')
        opt_ecx.par_read = 0;  % when set to 1 read parameters from saved set, but do NOT make CIs
        % this makes sure that calc_epx does not re-read the parameters every time
    end
    if type_conf < 1 || par_read == 1 % then we don't want to make CIs
        type_conf = 0; % just do analysis without CI
        rnd       = []; % make sample empty
    end
else 
    rnd = [];
end

%% Run through profiles
% This calls calc_epx_window_sub every time, for each profile. This is a
% function within this function.
% 
% Note: the parallel toolbox is used in calc_epx_window_sub to run calc_epx
% for CIs. Therefore, we cannot use parfor here, to run through profiles.
% In calc_epx_window_sub, the globals are still used. Advantage is that
% this code allows to batch-calculate CIs as well, without further trouble.

diary(glo.diary) % collect output in the diary "results.out"

opt_ecx.batch_epx = 1; % when set to 1 use batch mode (no output to screen)
% This makes sure that calc_epx does not produce output to screen.

COLL = [];
% If the profile is not located in the working directory, return the partial path
if ~strcmp(filepath(1:end-1),pwd)
    a = strfind(filepath,filesep);
    if length(a) > 2
        disp(['profiles from directory: ',filepath(a(end-2):end)])
    else
        disp(['profiles from directory: ',filepath])
    end
end

if type_conf < 1
    disp(['Calculating EPx without CIs, using moving time windows, for ',num2str(length(filename)),' exposure profiles.'])
else
    disp(['Calculating EPx with CIs, using moving time windows, for ',num2str(length(filename)),' exposure profiles.'])
end
disp('After finishing a profile, the results will be immediately shown on screen (and in the log-file).')
disp('After finishing all profiles, results are shown again as a table.')
disp('')

f = waitbar(0,'Calculating moving time window in batch mode.','Name','calc_epx_window_batch.m');

remX_coll = []; % collect the results for removed traits

for i = 1:length(filename)
    
    disp(['Calculating profile from file: ',filename{i},' (',num2str(i),' of ',num2str(length(filename)),')'])
        
    % <MinColl> collects, for each state, the MF, minimum of the trait relative
    % to the control (thus largest effect), and time at which it occurs.
    % <MinCI> collects the CI for the minimum trait level. Structure is
    % MinColl{i_MF}(i_X,:). Output of ind_traits is needed to know which state
    % is meant with i_X.
    [MinColl,MinCI,ind_traits,remX_trt] = calc_epx_window_sub(par_plot,[filepath,filename{i}],Twin,opt_ecx,opt_conf,rnd,glo);
    remX_coll = cat(1,remX_coll,remX_trt(:)); % this collects the position in the state variable list

    if ~isempty(ind_traits)
        for j = 1:length(Feff) % run through effect levels
            if ~iscell(MinColl) && isnan(MinColl)
                res_tmp = NaN; % catch cases where no trait has any effect in any window
                % this may not be needed, as ind_traits would be empty
            else
                res_tmp = MinColl{j};
            end
            if type_conf > 0
                ci_tmp  = MinCI{j};
            else
                ci_tmp  = [];
            end
            for k = 1:length(ind_traits) % run through endpoints
                if isempty(ci_tmp)
                    COLL = cat(1,COLL,[i res_tmp(k,:) ind_traits(k)]);
                else
                    COLL = cat(1,COLL,[i res_tmp(k,:) ind_traits(k) ci_tmp(k,:)]);
                end
            end
        end
        disp_epx(COLL(COLL(:,1)==i,:),filename(i),loc_all,mf_crit,calc_int,opt_ecx,type_conf,[],Twin,remX_trt)
    else
        disp(['  For this profile, no EPx (lower than max ',num2str(max_mf),') was found.'])
    end
    waitbar(i/length(filename),f) % update the waiting bar
end

close(f) % close the waiting bar

%% Display summary results on screen

if ~isempty(COLL)
    disp('Final table with results from all profiles.')
    disp_epx(COLL,filename,loc_all,mf_crit,calc_int,opt_ecx,type_conf,warn_init,Twin,remX_coll)
else
    disp('There are no results to display. This probably means that there was insufficient effect at the highest allowed multiplication factor, for all profiles.')
    disp(' ')
end

disp(['Time required: ' secs2hms(toc)])

diary off

function disp_epx(COLL,filename,loc_all,mf_crit,calc_int,opt_ecx,type_conf,warn_init,Twin,remX_trt)

rob_win   = opt_ecx.rob_win;   % set to 1 to use robust EPx calculation for moving time windows, rather than with fzero
rob_rng   = opt_ecx.rob_rng;   % range within which robust EPx is calculated, and number of points
max_mf    = opt_ecx.max_mf;    % maximum allowed multiplication factor
% disp(' ')

locS  = loc_all(1);
locL  = loc_all(2);
locR  = loc_all(3);
loc_a = loc_all(4);

Feff = unique(COLL(:,2));
warn = [0 0];

if ~isempty(warn_init) % if we only have 1 filename, we still want to see the end result with this header info
    disp('=================================================================================')
    disp('Results from batch calculations with moving time window')
    disp(['  Length of window  : ',num2str(Twin)])
    disp(['  stepsize window   : ',num2str(opt_ecx.Tstep)])
    disp(['  max MF calculated : ',num2str(max_mf)])

    switch type_conf
        case 0
            disp('  Lowest EPx, no confidence intervals');
        case 1
            disp('  Lowest EPx with CIs: Bayesian 95% credible interval');
        case 2
            disp('  Lowest EPx with CIs: 95% pred. likelihood, shooting method');
        case 3
            disp('  Lowest EPx with CIs: 95% pred. likelihood, parspace explorer');
    end
    if rob_win == 1
        disp(['  Robust EPx calculation with range between: ',num2str(rob_rng(1)),'-',num2str(rob_rng(end)),' (n=',num2str(length(rob_rng)),')'])
    end
    disp(' ')
end

for j = 1:length(Feff) % run through effect levels
    
    if length(filename) > 1
        disp('=================================================================================')
        disp(['Filename of profile       EP',num2str(Feff(j)*100)])
        disp('=================================================================================')
    else
        disp(['=========== EP',num2str(Feff(j)*100),' ========================================='])
        fprintf('%-30s \n',filename{1})
    end

    COLL_tmp = COLL(COLL(:,2)==Feff(j),:);
    
    for i = 1:length(filename)
        if length(filename) > 1
            fprintf('%-30s ',filename{i})
            ind_i = find(COLL_tmp(:,1)==i);
        else
            ind_i = (1:size(COLL_tmp,1))';
        end
        
        if isempty(ind_i) % for some profiles, there may be no effects, so ind_i is empty
            fprintf('no results to display (insufficient effects) \n')
        end

        for k = 1:length(ind_i) % run through endpoints
            if k > 1
                if length(filename) > 1
                    fprintf('%-30s ',' ')
%                 else
%                     fprintf('%-2s ',' ')
                end
            end

            switch COLL_tmp(ind_i(k),5)
                case -1
                    if calc_int == 1
                        fprintf('intr. rate ')
                    else
                        fprintf('surv-corr repro ')
                    end
                case locS
                    fprintf('%-7s ','surv')
                case locL
                    fprintf('%-7s ','length')
                case locR
                    fprintf('%-7s ','repro')
                case loc_a
                    fprintf('%-7s ','hlthy')
            end
            fprintf('%7.2f',COLL_tmp(ind_i(k),3))
            if length(COLL_tmp(ind_i(k),:)) > 5
                fprintf(' (%7.2f-%7.2f) ',COLL_tmp(ind_i(k),6:7))
            end
            fprintf(' t=%4.0f ',COLL_tmp(ind_i(k),4))
            if COLL_tmp(ind_i(k),3) < mf_crit
                fprintf('* ')
                warn(1) = 1;
            end
            if rob_win == 1
                if COLL_tmp(ind_i(k),3) == rob_rng(1) || COLL_tmp(ind_i(k),3) == rob_rng(end)
                    fprintf('#')
                    warn(2) = 1;
                end
            end
%             if length(filename) > 1
                fprintf('\n')
%             end

        end
%         if length(filename) == 1
%             fprintf('\n')
%         end
    end
end

disp('=================================================================================')
if warn(1) == 1
    disp(['* these EPx are below a critical value of ',num2str(mf_crit),' and may warrant closer look.'])
end
if warn(2) == 1
    disp('# the EPx is beyond the range that is searched; edge of range is reported.')
end
if ~isempty(warn_init)
    if warn_init(1) == 1
        disp('You are using separate parameters per data set (with glo.names_sep); make sure you use the right parameter set.')
    end
    if warn_init(2) == 1
        disp('For this combination of pMoA and feedbacks, robust calculation is advised to make certain that no lower EPx exist.')
    end
end
if ~isempty(remX_trt) % it is also possible to spell out the traits by comparing to locL, etc.
    disp('Some traits have been removed from the output as they produced no (or not enough) effect.')
end
disp(' ')

% if length(filename) > 1
%     % display again for easy copying and plotting TEMP TEMP
%     for j = 1:length(Feff) % run through effect levels
%         COLL_tmp = COLL(COLL(:,2)==Feff(j),:);
%         for i = 1:length(filename)
%             ind_i = find(COLL_tmp(:,1)==i);
%             for k = 1:length(ind_traits) % run through endpoints
%                 fprintf('%4.2f',COLL_tmp(ind_i(k),3))
%                 fprintf(' %4.0f ',COLL_tmp(ind_i(k),4))
%             end
%             for k = 1:length(ind_traits) % run through endpoints
%                 if length(COLL_tmp(ind_i(k),:)) > 5
%                     fprintf(' %4.2f %4.2f ',COLL_tmp(ind_i(k),6:7))
%                 end
%             end
%             fprintf('\n')
%         end
%     end
% end
% =========================================================================

function [MinColl,MinCI,ind_traits,remX_trt] = calc_epx_window_sub(par_plot,fname_prof,Twin,opt_ecx,opt_conf,rnd,glo)

% Syntax: [MinColl,MinCI,ind_traits] = calc_epx_window_sub(par_plot,fname_prof,Twin,opt_ecx,opt_conf,rnd,glo)
%          
% Calculate EPx/LPx due to an exposure profile, with moving time window.
% This function should work with every TKTD model you throw at it, as long
% as there is at least one state variable indicated with one of the
% dedicated traits (whose position in the state vector is given by:
% <glo.locS>, <glo.locL> or <glo.locR>). The actual calculations of the EPx
% (with CI) is performed by calc_epx. In this function, the CI is only
% calculated on the lowest EPx for each trait. Note: this is a simplified
% version of calc_epx_window.
% 
% <par_plot>   parameter structure for the best-fit curve; if left empty the
%              structure from the saved sample is used
% <fname_prof> filename for the file containing the exposure profile
% <Twin>       length of time window (days)
% <opt_ecx>    options structure for ECx and EPx calculations
% <opt_conf>   options structure for making confidence intervals
% <rnd>        the sample with parameter sets as input
% <glo>        the global structure glo, as input rather than global
% 
% <MinColl> collects, for each state and each effect level, the minimum
% EPx, and time at which it occurs. <MinCI> collects the CI for the minimum
% trait level. Structure is MinColl{i_eff}(i_X,:). Output of ind_traits is
% needed to know which state is meant with i_X.
% 
% Author     : Tjalling Jager 
% Date       : May 2023
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2022, Tjalling Jager.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE-MIT.txt file in the root directory of BYOM. 

X_excl    = opt_ecx.statsup;   % states to suppress from the calculations (e.g., locS)
start_neg = opt_ecx.start_neg; % set to 1 to start the moving window at minus window width
Feff      = opt_ecx.Feff;      % effect level (>0 en <1), x/100 in LCx (also used here for ECx)
prune_win = opt_ecx.prune_win; % set to 1 to prune the windows to keep the interesting ones
calc_int  = opt_ecx.calc_int;  % integrate survival and repro into 1) RGR, or 2) survival-corrected repro (experimental!)
Tstep     = opt_ecx.Tstep;     % stepsize or resolution of the time window (default 1 day)
saveall   = opt_ecx.saveall;   % set to 1 to save all output (EPx for each element of the sample) in separate folder
batch_plt = opt_ecx.batch_plt; % when set to 1 create plots when using calc_epx_window_batch

if isempty(opt_conf)
    type_conf = 0; % then we don't need CIs
else
    type_conf = opt_conf.type; % use values from slice sampler (1), likelihood region(2) to make intervals
    type_conf = max(0,type_conf); % if someone uses -1, set it to zero
end

% see which states are there
locS = [];
locL = [];
locR = [];
if isfield(glo,'locS') && ~ismember(glo.locS,X_excl) % then we have a state of survival
    locS = glo.locS; % collect the location
end
if isfield(glo,'locL') && ~ismember(glo.locL,X_excl) % then we have a state of body length
    locL = glo.locL; % collect the location
end
if isfield(glo,'locR') && ~ismember(glo.locR,X_excl) % then we have a state of reproduction
    locR = glo.locR; % collect the location
end
% And also add the states for the GUTS immobility package. For now, active
% only, since for that trait, it is easy to calculate ECx relative to the
% control (for death and immobile, the control is zero). There is a way to
% calculate EPx for death, but that would require summing active and
% immobile animals before calculating the effect (or take 1-death), so that
% is a bit more work.
loc_a = [];
if isfield(glo,'loc_a') % then we have a state of active
    loc_a = glo.loc_a; % collect the locations
end

ind_traits = [locS locL locR loc_a]; % indices for the traits we want from Xout
% We seem to have no interest for damage ...

% Note: removed code from here that loads the rnd, and placed it in main
% function. This function thus always gets par_plot and rnd as input.

% Load exposure profile from file to use with linear interpolation, and
% prepare it for moving-time window calculations
Cw = load(fname_prof);
if isempty(Twin) % if Twin is left empty, we just calculate the entire profile
    Trange = 0;
    Twin   = Cw(end,1); % last time point in profile
    indT   = [1 length(Cw)]; % indices for the start and stop point in Cw
else
    % Prepare exposure profile for moving-time window calculations
    [Cw,indT,Trange] = prep_exp_profile(Cw,Twin,Tstep,start_neg);
end

% figure
% plot(Cw(:,1),Cw(:,2),'k-')

% =================== TEST ================================================
if prune_win == 1
    Trange_sel = prune_windows(Trange,Cw,Twin); % prune the range to the interesting windows
end
% =================== TEST ================================================

opt_ecx.batch_epx = 1; % set to batch mode so calc_epx provides no output
opt_ecx.par_read  = 0; % when set to 1 read parameters from saved set, but do NOT make CIs
opt_conf.type     = type_conf;
% NOTE: I change the options, but not batch_epx etc.! That way, this
% funtion keeps the setting it was called with, and only calc_epx will go
% into batch mode.

N_traits = length(ind_traits);
if calc_int > 0 % if we integrate survival and repro ...
    N_traits   = 1; % we only have 1 trait left
    ind_traits = -1;
end

% Xcoll will collect all EPx output, and is initialised with NaNs
Xcoll = cell(1,length(Feff));
for j = 1:length(Feff)
    Xcoll{j} = nan(length(Trange),N_traits);
end

% Note: the actual calculations of the EPx (with CIs) are performed by
% calc_epx. Therefore, we cannot use parfor again here. It would be nice to
% use parfor for the different windows. However, this runs into problems
% with the scenarios, as make_scen is used in calc_epx, which stores the
% scenario in a global ... and within parallel processing, we cannot change
% globals. This requires a much more thought.

traitcoll = [];
for i_T = 1:length(Trange) % run through all time points (start of window)
    
    if prune_win == 1 && Trange_sel(i_T) == 0 % then this window is pruned
        continue % so move to next window!
    end
    
    T      = [Trange(i_T);Trange(i_T)+Twin]; % time window of length Twin
    Cw_tmp = Cw(indT(i_T,1):indT(i_T,2),:);  % extract only the profile that covers the time window
    if Cw_tmp(end,1)-Cw_tmp(1,1) ~= Twin     % to make absolutely sure ...
        error('Something went wrong: the exposure window length does not match the required length.')
    end
        
    Cw_tmp(:,1) = Cw_tmp(:,1)-T(1); % make time vector for the short profile start at zero again

    if ~all(Cw_tmp(:,2)==0) % if there is no exposure, there is no effect
        [EPx,~,~,ind_traits_tmp] = calc_epx(par_plot,Cw_tmp,[],opt_ecx,[],[]);
        % NOTE: not important to specify the time window with Twin here, as
        % Cw_tmp is modified to exactly match all time windows!
        traitcoll = unique([traitcoll ind_traits_tmp]); % collect traits that have meaningful output, at some windows
        for i_trt = 1:length(ind_traits_tmp)
            [~,ind_trt] = ismember(ind_traits_tmp(i_trt),ind_traits);
            for j_eff = 1:length(Feff)
                Xcoll{j_eff}(i_T,ind_trt) = EPx{j_eff}(i_trt);
            end
        end
    end
    
end

remX_trt = [];
if calc_int == 0
    [X_excl_add,remX] = setdiff(ind_traits,traitcoll); % traits that are not returned by calc_epx since effect is too small
    remX_trt = ind_traits(remX); % as extra output, use the position in the state variable list
    ind_traits(remX)  = []; % remove that trait from the trait list
    for j_eff = 1:length(Feff)
        Xcoll{j_eff}(:,remX) = []; % remove that trait from the collected values
    end
    N_traits = length(ind_traits); % redefine since traits may be removed
else 
    X_excl_add = [];
end

if isempty(traitcoll)
    MinColl    = NaN;
    MinCI      = NaN;
    ind_traits = [];
    return
end

%% Find where the lowest EPx occurs and how small it is

MinColl = cell(1,length(Feff));

for i_X = 1:N_traits % run through traits
    for i_eff = 1:length(Feff) % run through standard multiplication factors
        % find lowest EPx, and time window where it occurs, for each
        % state at each effect level.
        [val_min,ind_min] = min(Xcoll{i_eff}(:,i_X));
        
        MinColl{i_eff}(i_X,:) = [Feff(i_eff) val_min Trange(ind_min)]; % collect effect level, lowest EPx and time at which it occurs
    end
end

%% Now repeat calling calc_epx for the CIs in selected windows
% The CI calculation is very slow, so we only do it here for the specific
% windows where we found the lowest EPx (for each x and each trait).

if type_conf > 0

    MinCI = cell(1,length(Feff)); 

    opt_ecx.statsup = [X_excl(:);X_excl_add]; % states to suppress from the calculations
    % Need to redefine, since some may have produced no interesting output,
    % and were not returned by calc_EPx. So there is no need to ask to
    % calculate them again.

    for i_X = 1:N_traits % run through traits
        for i_eff = 1:length(Feff) % run through standard multiplication factors
            Tstrt = MinColl{i_eff}(i_X,3); % read start of crucial window for this effect level and trait
            T     = [Tstrt;Tstrt+Twin];    % time window of length Twin
            
            i_T    = find(Trange==Tstrt);           % index to start time in Cw
            Cw_tmp = Cw(indT(i_T,1):indT(i_T,2),:); % extract only the profile that covers the time window
            if Cw_tmp(end,1)-Cw_tmp(1,1) ~= Twin    % to make absolutely sure ...
                error('Something went wrong: the exposure window length does not match the required length.')
            end
            
            Cw_tmp(:,1) = Cw_tmp(:,1)-T(1); % make time vector for the short profile start at zero again
            if saveall == 0
                [~,EPx_lo,EPx_hi,ind_traits_tmp] = calc_epx(par_plot,Cw_tmp,[],opt_ecx,opt_conf,[],rnd);
            else
                [~,EPx_lo,EPx_hi,ind_traits_tmp] = calc_epx(par_plot,Cw_tmp,[],opt_ecx,opt_conf,[],rnd,Tstrt);
            end
            % NOTE: not important to specify the time window with Twin
            % here, as Cw_tmp is modified to exactly match all time windows!
            [~,ind_trt] = ismember(ind_traits(i_X),ind_traits_tmp); % where is our current trait in the output of calc_epx?
            
            MinCI{i_eff}(i_X,:) = [EPx_lo{i_eff}(ind_trt) EPx_hi{i_eff}(ind_trt)]; % collect that CI
        end
    end
else
    MinCI = NaN; % make sure output is defined
end

if batch_plt == 0 % when NOT making plots, we can return
    return
end

%% Plot the results
% This makes a multiplot: the different endpoints (traits) are shown in
% rows (EPx plotted), and the different effect levels in columns. Note that
% the x-axis is the time at which the time window starts; the window is
% shifted across the entire profile, until it extends beyond the profile.

showprof  = opt_ecx.showprof;  % set to 1 to show exposure profile in plots at top row
rob_win   = opt_ecx.rob_win;   % set to 1 to use robust EPx calculation for moving time windows, rather than with fzero
rob_rng   = opt_ecx.rob_rng;   % range within which robust EPx is calculated, and number of points
lim_yax   = opt_ecx.lim_yax;   % limit for y-axiswith EPx values (calc_epx_window)
notitle   = opt_ecx.notitle;   % set to 1 to suppress titles above plots

n = length(Feff);
m = length(ind_traits);
if showprof > 0
    m = m + 1;
end

[figh,ft] = make_fig(m,n,2); % create a figure window of correct size

if showprof > 0

    for i_eff = 1:length(Feff) % run through effect levels
        h_pl = subplot(m,n,i_eff); % make a sub-plot
        hold on
        h_ax = gca;
        set(h_ax,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting
        if m>1 && n>1 % only shrink white space when there are more than 1 rows and columns
            p = get(h_pl,'position');
            p([3 4]) = p([3 4])*1.1; % add 10 percent to width and height
            set(h_pl, 'position', p);
        end

        if i_eff == 1
            h_tmp = ylabel('Exposure conc.'); set(h_tmp,ft.name,ft.label)
        else
            set(h_ax,'YTickLabel',[]); % remove tick labels on y-axis
        end
        set(h_ax,'XTickLabel',[]); % remove tick labels on x-axis
        % h_tmp = xlabel(glo.xlab); set(h_tmp,ft.name,ft.label)
        % 

        mod_t = Cw(:,1);
        mod_plt(:,1) = Cw(:,2);
        mod_plt(:,2) = zeros(size(mod_plt,1),1);
        plot(mod_t,mod_plt,'k-')

        % Little trick to fill the area between the two curves, to
        % obtain a coloured confidence interval as a band.
        t2  = [mod_t;flipud(mod_t)]; % make a new time vector that is old one, plus flipped one
        Xin = [mod_plt(:,1);flipud(mod_plt(:,2))]; % do the same for the plot line, hi and lo
        fill(t2,Xin,'b','LineStyle','none','FaceAlpha',1) % and fill this object

        a = ylim;
        plot([0 0],[0 a(2)],'k:') % plot vertical line for start of profile

        xlim([Trange(1) Trange(end)]) % limit x-axis
        title(['Feff = ',num2str(Feff(i_eff))])

    end
end

for i_X = 1:N_traits % run through traits
    
    for i_eff = 1:length(Feff) % run through effect levels
        
        if showprof == 0
            h_pl = subplot(m,n,(i_X-1)*length(Feff)+i_eff); % make a sub-plot
        else
            h_pl = subplot(m,n,(i_X)*length(Feff)+i_eff); % make a sub-plot
        end
        hold on
        h_ax = gca;
        set(h_ax,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting
        if m>1 && n>1 % only shrink white space when there are more than 1 rows and columns
            p = get(h_pl,'position');
            p([3 4]) = p([3 4])*1.1; % add 10 percent to width and height
            set(h_pl, 'position', p);
        end
        
        % create axis labels
        if i_X == length(ind_traits) % only x-label at last row
            xlab = ['start ' glo.xlab]; % label
            h_tmp = xlabel(xlab); set(h_tmp,ft.name,ft.label)
        else
            set(h_ax,'XTickLabel',[]); % remove tick labels on x-axis
        end
        if i_eff == 1 % only y-label in first column
            switch calc_int
                case 0
                    ylab = glo.ylab{ind_traits(i_X)};
                case 1
                    ylab = 'relative intrinsic rate';
                case 2
                    ylab = 'survival-corrected repro';
            end
            
            if rob_win == 1 % when using the robust routine
                ylab = ['rEPx ' ylab]; %#ok<AGROW> label
            else
                ylab = ['EPx ' ylab]; %#ok<AGROW> label
            end
            h_tmp = ylabel(ylab); set(h_tmp,ft.name,ft.label)
        else
            set(h_ax,'YTickLabel',[]); % remove tick labels on y-axis
        end
        if i_X == 1 && showprof == 0 % only title in first row
            title(['Feff = ',num2str(Feff(i_eff))])
        end
        
        if length(Trange) == 1 % have to come up with something special for when there is only 1 time point
            
            % also plot critical lines for safety margins of 10 and 100
            plot([-1 1],[10 10],'k--')
            plot([-1 1],[100 100],'k--')
            plot([-1 1],[1 1],'k--')
                 
            if type_conf > 0 % then we have CIs to plot on the lowest EPx only!
                errorbar(0,MinColl{i_eff}(i_X,2),MinColl{i_eff}(i_X,2)-MinCI{i_eff}(i_X,1),MinCI{i_eff}(i_X,2)-MinColl{i_eff}(i_X,2),'ko-','MarkerFaceColor','r','LineWidth',1)
            else
                plot(0,MinColl{i_eff}(i_X,2),'ko','MarkerFaceColor','r')
            end
            
            ylim([0.1 lim_yax]) % limit y-axis
            xlim([-1 1]) % limit x-axis
            set(gca, 'YScale', 'log')

        else

            plot(Trange,Xcoll{i_eff}(:,i_X),'k-','LineWidth',1) % plot effects relative to control
            if rob_win == 1 % when using the robust routine, there are hard min-max bounds
                ylim([rob_rng(1) rob_rng(end)]) % limit y-axis
                if rob_rng(1)>0
                    set(gca, 'YScale', 'log')
                end
            else
                ylim([0.1 lim_yax]) % limit y-axis
                set(gca, 'YScale', 'log')
            end
            xlim([Trange(1) Trange(end)]) % limit x-axis

            % also plot critical lines for safety margins of 1, 10 and 100
            plot([Trange(1) Trange(end)],[100 100],'k--')
            plot([Trange(1) Trange(end)],[10 10],'k--')
            plot([Trange(1) Trange(end)],[1 1],'k--')
            plot([0 0],[0.1 lim_yax],'k:') % and for start of profile

            % plot the lowest EPx (with its CI) on top of everything else
            if type_conf > 0 % then we have CIs to plot on the lowest EPx only!
                errorbar(MinColl{i_eff}(i_X,3),MinColl{i_eff}(i_X,2),MinColl{i_eff}(i_X,2)-MinCI{i_eff}(i_X,1),MinCI{i_eff}(i_X,2)-MinColl{i_eff}(i_X,2),'ko-','MarkerFaceColor','r','LineWidth',1)
            else
                plot(MinColl{i_eff}(i_X,3),MinColl{i_eff}(i_X,2),'ko','MarkerFaceColor','r')
            end

        end
        
    end
end

if notitle == 0
    switch type_conf
        case 0
            tit = 'No confidence intervals';
        case 1
            tit = 'CIs: Bayesian 95% credible interval';
        case 2
            tit = 'CIs: 95% pred. likelihood, shooting method';
        case 3
            tit = 'CIs: 95% pred. likelihood, parspace explorer';
    end
    axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping','off');
    text(0.5, 1,tit,'HorizontalAlignment','center','VerticalAlignment', 'top');
end

drawnow

if glo.saveplt > 0 % if we want to save the plot
    a = strfind(fname_prof,filesep);
    savenm = ['epx_window_plot_batch_',fname_prof(a(end)+1:end-4)];%
    save_plot(figh,savenm);
end
