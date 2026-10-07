function [MinColl,MinCI,ind_traits] = calc_epx_window(par_plot,fname_prof,Twin,opt_ecx,opt_conf)

% Syntax: [MinColl,MinCI,ind_traits] = calc_epx_window(par_plot,fname_prof,Twin,opt_ecx,opt_conf)
% 
% Calculate EPx/LPx due to an exposure profile, with moving time window.
% This function should work with every TKTD model you throw at it, as long
% as there is at least one state variable indicated with one of the
% dedicated traits (whose position in the state vector is given by:
% <glo.locS>, <glo.locL> or <glo.locR>). The actual calculations of the EPx
% (with CI) is performed by calc_epx.
% 
% <par_plot>   parameter structure for the best-fit curve; if left empty the
%              structure from the saved sample is used
% <fname_prof> filename for the file containing the exposure profile
% <Twin>       length of time window (days)
% <opt_ecx>    options structure for ECx and EPx calculations
% <opt_conf>   options structure for making confidence intervals
% 
% <MinColl> collects, for each state and each effect level, the minimum
% EPx, and time at which it occurs. <MinCI> collects the CI for the minimum
% trait level. Structure is MinColl{i_eff}(i_X,:). Output of ind_traits is
% needed to know which state is meant with i_X.
% 
% Author     : Tjalling Jager 
% Date       : July 2023
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2024, Tjalling Jager.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE-MIT.txt file in the root directory of BYOM. 

global glo

filenm    = glo.basenm;

X_excl    = opt_ecx.statsup;   % states to suppress from the calculations (e.g., locS)
par_read  = opt_ecx.par_read;  % when set to 1 read parameters from saved set, but do NOT make CIs
batch_epx = opt_ecx.batch_epx; % when set to 1 use batch mode (no output to screen)
notitle   = opt_ecx.notitle;   % set to 1 to suppress titles above plots
start_neg = opt_ecx.start_neg; % set to 1 to start the moving window at minus window width
Feff      = opt_ecx.Feff;      % effect level (>0 en <1), x/100 in LCx (also used here for ECx)
rob_win   = opt_ecx.rob_win;   % set to 1 to use robust EPx calculation for moving time windows, rather than with fzero
rob_rng   = opt_ecx.rob_rng;   % range within which robust EPx is calculated, and number of points
prune_win = opt_ecx.prune_win; % set to 1 to prune the windows to keep the interesting ones
calc_int  = opt_ecx.calc_int;  % integrate survival and repro into 1) RGR, or 2) survival-corrected repro (experimental!)
Tstep     = opt_ecx.Tstep;     % stepsize or resolution of the time window (default 1 day)
id_sel    = opt_ecx.id_sel;    % scenario to use from X0mat, scenario identifier, flag for ECx to use scenarios rather than concentrations
saveall   = opt_ecx.saveall;   % set to 1 to save all output (EPx for each element of the sample) in separate folder
showprof  = opt_ecx.showprof;  % set to 1 to show exposure profile in plots at top row
lim_yax   = opt_ecx.lim_yax;   % limit for y-axiswith EPx values (calc_epx_window)
max_mf    = opt_ecx.max_mf;   % maximum allowed multiplication factor

if isempty(opt_conf)
    type_conf = 0; % then we don't need CIs
else
    type_conf = opt_conf.type; % use values from slice sampler (1), likelihood region(2) to make intervals
    type_conf = max(0,type_conf); % if someone uses -1, set it to zero
end

if batch_epx == 0 % don't show warnings if we're in batch mode
    diary(glo.diary) % collect output in the diary "results.out"
    warning('off','backtrace')
    if ~isempty(glo.names_sep)
        warning('You are using separate parameters per data set (with glo.names_sep)')
        warning(['Parameters for data set ',num2str(1+floor(id_sel(2)/100)),' will be used for EPx'])
        disp(' ')
    end
    if rob_win == 0 && isfield(glo,'moa') && isfield(glo,'feedb')
        if glo.feedb(2) == 1 && any(glo.moa(1:3)==1)
            warning('Calculating standard EPx windows is NOT recommended for this combination of pMoA and feedbacks!')
            warning('There is a possibility of (limited) lower EPx values than those reported.')
            warning('It is advisable to use the robust settings with opt_ecx.rob_win = 1.')
            disp(' ')
        end
    end
    warning('on','backtrace')
    diary off
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

% If we need CIs, load the best parameter set and the random sample from file
if type_conf > 0 || isempty(par_plot) % also if par_plot is not provided
    [rnd,par] = load_rnd(opt_conf); % loading and selecting the sample is handled in a separate function
    if numel(rnd) <= 1 % that means that no (useful) sample was found
        type_conf = 0; % no need to produce an error, just do analysis without CI
    end
    if isempty(par_plot) % if no par structure was entered in this function ...
        par_plot = par; % simply use the one from the sample file
    end
    if type_conf < 1 || par_read == 1 % then we don't want to make CIs
        type_conf = 0;  % don't make CIs anymore (when triggered by par_read)
        rnd       = []; % make sample empty
    end
else
    rnd = [];
end

% Load exposure profile from file to use with linear interpolation, and
% prepare it for moving-time window calculations
Cw     = load(fname_prof);
Tend   = Cw(end,1); % last time point in profile

if isempty(Twin) % if Twin is left empty, we just calculate the entire profile
    Trange = 0;
    Twin   = Tend; % last time point in profile
    indT   = [1 length(Cw)]; % indices for the start and stop point in Cw
else
    % Prepare exposure profile for moving-time window calculations
    [Cw,indT,Trange] = prep_exp_profile(Cw,Twin,Tstep,start_neg);
end
Tstrt = Cw(1,1);

% figure
% plot(Cw(:,1),Cw(:,2),'k-')

% =================== TEST ================================================
if prune_win == 1
    Trange_sel = prune_windows(Trange,Cw,Twin); % prune the range to the interesting windows
end
% =================== TEST ================================================

if batch_epx == 0 % don't show waitbar if we're in batch mode
    if rob_win == 0
        if type_conf > 0 % if we make CIs ...
            f = waitbar(0,'Calculating EPx with CIs for moving time window. Please wait.','Name','calc_epx_window.m');
        else
            f = waitbar(0,'Calculating EPx for moving time window. Please wait.','Name','calc_epx_window.m');
        end
    else
        if type_conf > 0 % if we make CIs ...
            f = waitbar(0,'Calculating robust EPx with CIs for moving time window. Please wait.','Name','calc_epx_window.m');
        else
            f = waitbar(0,'Calculating robust EPx for moving time window. Please wait.','Name','calc_epx_window.m');
        end
    end
end
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

Xcoll = cell(1,length(Feff));
Xlo   = cell(1,length(Feff));
Xhi   = cell(1,length(Feff));
% Xcoll will collect all EPx output, and is initialised with NaNs
for j = 1:length(Feff)
    Xcoll{j} = nan(length(Trange),N_traits);
    if type_conf > 0
        Xlo{j} = nan(length(Trange),N_traits);
        Xhi{j} = nan(length(Trange),N_traits);
    end
end

% Note: the actual calculations of the EPx (with CIs) are performed by
% calc_epx.

traitcoll = [];
for i_T = 1:length(Trange) % run through all time points (start of window)
    
    if batch_epx == 0
        waitbar(i_T/length(Trange),f) % update the waiting bar
    end
    
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
        if saveall == 0
            [EPx,EPx_lo,EPx_hi,ind_traits_tmp] = calc_epx(par_plot,Cw_tmp,[],opt_ecx,opt_conf,[],rnd);
        else % calling calc_epx with the start time for the window triggers saving of results!
            [EPx,EPx_lo,EPx_hi,ind_traits_tmp] = calc_epx(par_plot,Cw_tmp,[],opt_ecx,opt_conf,[],rnd,Trange(i_T));
        end
        % NOTE: not important to specify the time window with Twin here, as
        % Cw_tmp is modified to exactly match all time windows!
        traitcoll = unique([traitcoll ind_traits_tmp]); % collect traits that have meaningful output, at some windows
        for i_trt = 1:length(ind_traits_tmp)
            [~,ind_trt] = ismember(ind_traits_tmp(i_trt),ind_traits);
            for j_eff = 1:length(Feff)
                Xcoll{j_eff}(i_T,ind_trt) = EPx{j_eff}(i_trt);
                if type_conf > 0
                    Xlo{j_eff}(i_T,ind_trt) = EPx_lo{j_eff}(i_trt);
                    Xhi{j_eff}(i_T,ind_trt) = EPx_hi{j_eff}(i_trt);
                end
            end
        end
    end
    
end

if batch_epx == 0
    close(f) % close the waiting bar
end

if calc_int == 0
    
    [~,remX] = setdiff(ind_traits,traitcoll); % traits that are not returned by calc_epx since effect is too small
    
    if batch_epx == 0 && ~isempty(remX)
        disp(' ')
        for i = 1:length(remX)
            switch ind_traits(remX(i))
                case locS
                    disp('For survival, there is insufficient range of effects to calculate all EPx.')
                case locL
                    disp('For body length, there is insufficient range of effects to calculate all EPx.')
                case locR
                    disp('For reproduction, there is insufficient range of effects to calculate all EPx.')
                case loc_a
                    disp('For active animals, there is insufficient range of effects to calculate all EPx.')
            end
        end
    end

    ind_traits(remX) = []; % remove that trait from the trait list
    for j_eff = 1:length(Feff)
        Xcoll{j_eff}(:,remX) = []; % remove that trait from the collected values
        if type_conf > 0
            Xlo{j_eff}(:,remX) = []; % remove that trait from the collected values
            Xhi{j_eff}(:,remX) = []; % remove that trait from the collected values
        end
    end
    N_traits = length(ind_traits); % redefine since traits may be removed
    
end

if isempty(traitcoll)
    MinColl    = NaN;
    MinCI      = NaN;
    ind_traits = [];
    disp(['Stopping here, no effects found for any trait at max MF: ',num2str(max_mf)])
    return
end

%% Find where the lowest EPx occurs and how small it is

MinColl = cell(1,length(Feff));
MinCI   = cell(1,length(Feff));
for i_X = 1:N_traits % run through traits
    for i_eff = 1:length(Feff) % run through standard multiplication factors
        % find lowest EPx, and time window where it occurs, for each
        % state at each effect level.
        [val_min,ind_min] = min(Xcoll{i_eff}(:,i_X));
        
        MinColl{i_eff}(i_X,:) = [Feff(i_eff) val_min Trange(ind_min)]; % collect effect level, lowest EPx and time at which it occurs
        if type_conf > 0
            MinCI{i_eff}(i_X,:)   = [Xlo{i_eff}(ind_min,i_X) Xhi{i_eff}(ind_min,i_X)]; % also collect CI on the lowest EPx
        end
    end
end

if batch_epx == 1  % in batch mode, we can stop here
    return 
end

%% Plot the results
% This makes a multiplot: the different endpoints (traits) are shown in
% rows (EPx plotted), and the different effect levels in columns. Note that
% the x-axis is the time at which the time window starts; the window is
% shifted across the entire profile, until it extends beyond the profile.

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
        % xlabel(glo.xlab,ft.name,ft.label)
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

        xlim([Tstrt Tend]) % limit x-axis
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
                        
            errorbar(0,Xcoll{i_eff}(:,i_X),Xcoll{i_eff}(:,i_X)-Xlo{i_eff}(:,i_X),Xhi{i_eff}(:,i_X)-Xcoll{i_eff}(:,i_X),'ko-','MarkerFaceColor','k','LineWidth',1)
            ylim([0.1 lim_yax]) % limit y-axis
            xlim([-1 1]) % limit x-axis
            set(gca, 'YScale', 'log')

        else
            
            if type_conf > 0 % then we have CIs to plot
                % Little trick to fill the area between the two curves, to
                % obtain a coloured confidence interval as a band.
                t2  = [Trange';flipud(Trange')]; % make a new time vector that is old one, plus flipped one
                Xin = [Xlo{i_eff}(:,i_X);flipud(Xhi{i_eff}(:,i_X))]; % do the same for the plot line, hi and lo
                fill(t2(~isnan(Xin)),Xin(~isnan(Xin)),'g','LineStyle','none','FaceAlpha',1) % and fill this object
            end
            
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
            xlim([Tstrt Tend]) % limit x-axis
            
            if type_conf > 0 % then we have CIs to plot
                plot(Trange,Xlo{i_eff}(:,i_X),'k:') % plot effects
                plot(Trange,Xhi{i_eff}(:,i_X),'k:') % plot effects
            end
            
            % also plot critical lines for safety margins of 1, 10 and 100
            plot([Tstrt Tend],[100 100],'k--')
            plot([Tstrt Tend],[10 10],'k--')
            plot([Tstrt Tend],[1 1],'k--')
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
    savenm = ['epx_window_plot_',filenm];%
    save_plot(figh,savenm);
end

%% Print results on screen

% We might be called with an entire path to a file, rather than just a
% filename. In that case, only display the filename with the lowest folder
% it is in.
a = strfind(fname_prof,filesep);
if length(a) > 2
    fname_prof = fname_prof(a(end-1):end);
end

diary(glo.diary) % collect output in the diary "results.out"
% disp(' ')
disp('Results from moving time window calculations')
disp(['  Exposure from file: ',fname_prof])
disp(['  Length of window  : ',num2str(Twin)])
disp(['  stepsize window   : ',num2str(Tstep)])
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
disp('================================================================================')
disp('  Lowest EPx (most sensitive time window) for various traits')
disp('  Note that multiple minima may occur; only first one is indicated')
disp('  Also note that most sensitive window may depend on the effect level')
if rob_win == 1
    disp(['  Robust EPx calculation with range between: ',num2str(rob_rng(1)),'-',num2str(rob_rng(end)),' (n=',num2str(length(rob_rng)),')'])
end
disp('================================================================================')
for i_eff = 1:length(Feff) % run through effect levels
    disp(['EP',num2str(100*Feff(i_eff))])
    for i_X = 1:N_traits
        if calc_int > 0
            if calc_int == 1
                fprintf('  Intrinisic rate: ')
            else
                fprintf('  Surv-corr repro: ')
            end
        else
            switch ind_traits(i_X)
                case locS
                    fprintf('  Survival    : ')
                case locL
                    fprintf('  Body length : ')
                case locR
                    fprintf('  Reproduction: ')
                case loc_a
                    fprintf('  Active      : ')
            end
        end
                
        % Assume that time window always starts on a full day (the %3.0f
        % only prints whole numbers). Here, I convert the values into
        % strings; this allows me to use < and > for the robust EPx, which
        % has a hard boundary.
        a1 = sprintf('%#.2f',MinColl{i_eff}(i_X,2));
        if type_conf > 0
            a2 = sprintf('%#.2f',MinCI{i_eff}(i_X,1));
            a3 = sprintf('%#.2f',MinCI{i_eff}(i_X,2));
        end
        if rob_win == 1
            if MinColl{i_eff}(i_X,2)==rob_rng(1) 
                a1 = sprintf('<%#.2f',MinColl{i_eff}(i_X,2));
            elseif MinColl{i_eff}(i_X,2)==rob_rng(end)
                a1 = sprintf('>%#.2f',MinColl{i_eff}(i_X,2));
            end
            if type_conf > 0
                if MinCI{i_eff}(i_X,1)==rob_rng(1)
                    a2 = sprintf('<%#.2f',MinCI{i_eff}(i_X,1));
                end
                if MinCI{i_eff}(i_X,2)==rob_rng(end)
                    a3 = sprintf('>%#.2f',MinCI{i_eff}(i_X,2));
                end
            end
        end
        fprintf('%10s at t = %3.1f ',a1,MinColl{i_eff}(i_X,3))
        if type_conf > 0
            fprintf('(%s - %s)',a2,a3)
        end
        fprintf('\n')
    end
end
disp('================================================================================')
if type_conf > 0
    disp('Note: the CIs on the EPx are the CIs for the EPx in the time window with the lowest best-estimate EPx.')
end
if ~isempty(remX)
    disp('Some traits have been removed from the output as they produced no (or not enough) effect.')
end

if type_conf > 0
    disp(' ')
    disp('================================================================================')

    for i_eff = 1:length(Feff) % run through effect levels
        disp(['Estimate with lowest minimal CI bound for EP',num2str(100*Feff(i_eff)),' across all windows'])
        for i_X = 1:N_traits
            if calc_int > 0
                if calc_int == 1
                    fprintf('  Intrinisic rate: ')
                else
                    fprintf('  Surv-corr repro: ')
                end
            else
                switch ind_traits(i_X)
                    case locS
                        fprintf('  Survival    : ')
                    case locL
                        fprintf('  Body length : ')
                    case locR
                        fprintf('  Reproduction: ')
                    case loc_a
                        fprintf('  Active      : ')
                end
            end
            [minEPx,indEPx] = min(Xlo{i_eff}(:,i_X),[],1,'omitnan');
            bestEPx = Xcoll{i_eff}(indEPx,i_X); % best EPx belonging to lowest CI
            maxEPx = Xhi{i_eff}(indEPx,i_X);

            % Assume that time window always starts on a full day (the %3.0f
            % only prints whole numbers). Here, I convert the values into
            % strings; this allows me to use < and > for the robust EPx, which
            % has a hard boundary.
            a1 = sprintf('%#.2f',bestEPx);
            a2 = sprintf('%#.2f',minEPx);
            a3 = sprintf('%#.2f',maxEPx);
            if rob_win == 1
                if bestEPx==rob_rng(1)
                    a1 = sprintf('<%#.2f',bestEPx);
                elseif bestEPx==rob_rng(end)
                    a1 = sprintf('>%#.2f',bestEPx);
                end
                if minEPx==rob_rng(1)
                    a2 = sprintf('<%#.2f',minEPx);
                end
                if maxEPx==rob_rng(end)
                    a3 = sprintf('>%#.2f',maxEPx);
                end
            end
            fprintf('%10s at t = %3.1f ',a1,Trange(indEPx))
            if type_conf > 0
                fprintf('(%s - %s)',a2,a3)
            end
            fprintf('\n')
        end
    end
    disp('================================================================================')
end

disp(['Time required: ' secs2hms(toc)])
diary off

% % minimum window for repro across effect levels
% MinRep = [];
% for i_eff = 1:length(Feff) % run through standard multiplication factors
%     for i_X = 1:length(ind_traits)
%         if ind_traits(i_X) == locR % then we have repro data
%             % collect MF, relative cumul. repro, and start of time window
%             MinRep = cat(1,MinRep,[Feff(i_eff) MinColl{i_eff}(i_X,2) MinColl{i_eff}(i_X,3)]);
%         end
%     end
% end
% % How to find the 'most sensitive time window'? The window with most effect
% % WILL depend on the value of the MF. Furthermore, for different elements
% % of the sample, the most-sensitive window will be somewhere else.
% % Furthermore, it will differ between endpoints.



