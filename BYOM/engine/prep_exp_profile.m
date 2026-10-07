function [Cw_out,ind_out,Trange] = prep_exp_profile(Cw,Twin,Tstep,start_neg)

% Usage: [Cw_out,ind_out,Trange] = prep_exp_profile(Cw,Twin,Tstep,start_neg)
% 
% Small helper function to prepare an exposure profile for specific- or
% moving-window calculation. This function is called from calc_epx,
% calc_epx_window and calc_epx_window_batch. It adds zeros at the end of
% the profile, makes sure that all start and stop time points are included
% in the profile (interpolated linearly when needed), and delivers a
% helpful matrix ind_out with the start-stop indices for each time window.
% That should save some calculation time, as finding start/stop points in
% the long matrix Cw is not done in a for loop for every time window
% anymore.
% 
% Author     : Tjalling Jager 
% Date       : May 2023
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2024, Tjalling Jager.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE-MIT.txt file in the root directory of BYOM. 

if length(Twin) == 1 % when called by calc_epx_window or calc_epx_window_batch

    Tend   = Cw(end,1); % last time point in profile (before we start adding points at the end)

    % Have the option to start with windows *before* the exposure profile.
    % This would be needed if the pulses come very early in the profile
    % (otherwise, they would hit juveniles only). This adds two one points,
    % one for the first window start, and one close to zero at the smallest
    % difference in the time vector of Cw.
    if start_neg == 1
        Cw = cat(1,[-1*Tstep*floor(Twin/Tstep) 0;-1*min(diff(Cw(:,1))) 0],Cw);
        % Note: I use Tstep*floor(Twin/Tstep) rather than Twin. This is for
        % cases where Twin is not a multiple of Tstep (which would probably
        % be rare). This ensures that there is always a window starting at
        % t=0.
    end

    % We may also want to include time windows that start almost at the end
    % of the profile, and thus that extend longer than the profile itself.
    % We can solve that by adding two time points after the profile that
    % are zero. The first point is at the minimum difference in the time
    % vector of Cw. The last point is probably not needed, but ensures that
    % interpolation always functions properly.
    Cw = cat(1,Cw,[Cw(end,1)+min(diff(Cw(:,1))) 0;Cw(end,1)+Twin 0]);

    Tstrt  = Cw(1,1);   % first time point in new profile
    Trange = Tstrt:Tstep:Tend; % start the time window every step, starting at start of the profile
    % Note: if Tend-Tstrt is not a multiple of Tstep, the last window will
    % start before the end of Cw. This is fine; the next window would start
    % after the profile, so be all zeros.

    % Modify the exposure profile to ensure that all start and end times, for
    % all windows, are included. Missing points are added with linear
    % interpolation.
    Tall   = unique([Trange,Trange+Twin]);
    Tint   = setdiff(Tall,Cw(:,1)); % find if start/end points are NOT in Cw
    Cw_new = interp1(Cw(:,1),Cw(:,2),Tint); % interpolate to the exact missing points in the profile
    Cw_out = [Cw;[Tint' Cw_new']]; % add the interpolated points to the exposure profile
    Cw_out = sortrows(Cw_out,1); % and sort on first column (time)
    
    [~,ind_1] = ismember(Trange,Cw_out(:,1));      % indices to start times in Cw
    [~,ind_2] = ismember(Trange+Twin,Cw_out(:,1)); % indices to stop times in Cw
    ind_out   = [ind_1' ind_2'];

else % when called by calc_epx

    % We may also want to include time windows that start almost at the end of
    % the profile, and thus that extend longer than the profile itself. We can
    % solve that by adding two time points after the profile that are zero. The
    % last point is probably not needed, but it does not hurt either.
    Cw = cat(1,Cw,[Cw(end,1)+min(diff(Cw(:,1))) 0;Cw(end,1)+diff(Twin) 0]);

    % Modify the exposure profile to ensure that all start and end times, for
    % all windows, are included. Missing points are added with linear
    % interpolation.
    Tint    = setdiff(Twin,Cw(:,1)); % find if start/end points are NOT in Cw
    Cw_new  = interp1(Cw(:,1),Cw(:,2),Tint); % interpolate to the exact points in the profile
    Cw_out  = [Cw;[Tint' Cw_new']]; % add the interpolated points to the exposure profile
    Cw_out  = sortrows(Cw_out,1); % and sort on first column (time)
    
    [~,ind_1]   = ismember(Twin(1),Cw_out(:,1)); % indices to start time in Cw
    [~,ind_2]   = ismember(Twin(2),Cw_out(:,1)); % indices to stop time in Cw
    Cw_out      = Cw_out(ind_1:ind_2,:); % only keep the range within Twin
    Cw_out(:,1) = Cw_out(:,1)-Twin(1);  % make time vector for the short profile start at zero again
    
    ind_out  = [];
    Trange   = [];
end