%% BYOM, combine_mat_files.m
%
% * Author: Tjalling Jager
% * Date: December 2021
% * Web support: <http://www.debtox.info/byom.html>
%
% *This script:* part of the package for mixture analyses with GUTS. This
% script asks the user to select two previously generated MAT files from
% the parameter-space explorer. These will be combined into parameter
% clouds representing the binary mixture under damage addition and
% independent action (new MAT files will be saved, that can later be used
% to create predictions and comparisons to mixture data). Note that this
% analysis assumes that there are NO interactions in the mixture.
% 
% Note: at this moment, the MAT files need to be in the same working
% directory as this script.
% 
% At this moment, the script can run through SD and IT sequentially, and
% only allows you to select the files with the correct indicator for the
% death mechanism.

%  Copyright (c) 2012-2021, Tjalling Jager.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE-MIT.txt file in the root directory of BYOM. 

clear, clear global % clear the workspace and globals

global glo2

pathdefine(0) % set path to the BYOM/engine directory (option 1 uses parallel toolbox)

% Note: the check for number of cores below is copied from prelim_checks.m.
% Since that function is skipped here, we need to have this code here as
% well. It is only used when calling pathdefine with the option 1, and when
% the parallel toolbox is installed.

% This checks whether you have the parallel toolbox installed. If it is
% installed, it calculates how many cores to use for parallel computation.
% Note that <parfor> loops will be executed as <for> loops without this
% toolbox. 
% 
% NOTE: the code below is copied from prelim_checks in engine_par. We need
% this since this script will be called without using prelim_checks.
n_cores = 0; % by default, we have zero cores
% this is to avoid errors when running this code without the toolbox
v = ver; % take the version information
if any(strcmp('Parallel Computing Toolbox',{v.Name}))
    % If a parallel pool has not yet been created, create one with n_cores-1 workers
    n_cores = feature('numcores');
    if n_cores < 7 % number of physical cores available
        n_cores = n_cores-1; % this leaves 1 core for other work
    else % if you have at least 8 cores, leave 2 cores free
        n_cores = n_cores-2; % this leaves 2 cores for other work
        % on my 10-core workstation, using 9 does not improve speed over
        % using 8 cores. Using 8 beats using 6 cores, though.
    end
    % Note: the pool is not automatically deleted after the analysis, but
    % closes after 30 min idle, by default. Otherwise, kill it manually:
      % poolobj = gcp('nocreate'); % If no pool, do not create new one.
      % delete(poolobj) % force parallel pool to end
end
glo2.n_cores = n_cores; % put the number of cores to use in glo2

SEL = [1 2]; % These are the death mechanisms that will be run automatically
% Note: option 3 (both SD and IT, so GUTS proper) is not implemented.

% Ask for filenames for the samples. Let's assume that the SD and IT files
% have the same name, apart from the 'sel*' part. NOTE: this only allows
% selection of files that start with "byom_" here! That is handy to avoid
% seeing the "ADD_" and "IND_" files that may be in the same folder.
filenmA = uigetfile(['byom_*sel',num2str(SEL(1)),'_PS.mat'],'Select saved project file for chemical A for comparison of parameter space','MultiSelect','off'); % use Matlab GUI to select MAT file
if ~iscell(filenmA) && numel(filenmA) == 1 && filenmA(1) == 0 % if cancel is pressed ...
    return % simply stop as there is no useful input to work with
end
filenmB = uigetfile(['byom_*sel',num2str(SEL(1)),'_PS.mat'],'Select saved project file for chemical B for comparison of parameter space','MultiSelect','off'); % use Matlab GUI to select MAT file
if ~iscell(filenmB) && numel(filenmB) == 1 && filenmB(1) == 0 % if cancel is pressed ...
    return % simply stop as there is no useful input to work with
end

for i = 1:length(SEL)
    if ~ismember(SEL(i),[1 2])
        error('Unknown death mechanism requested')
    end
       
    filename{1} = [filenmA(1:end-8),num2str(SEL(i)),'_PS.mat'];
    filename{2} = [filenmB(1:end-8),num2str(SEL(i)),'_PS.mat'];
    
    if exist(filename{1},'file') ~= 2
        error(['File <',filename{1},'> not found'])
    end
    if exist(filename{2},'file') ~= 2
        error(['File <',filename{2},'> not found'])
    end
    if strcmp(filename{1},filename{2})
        error('Chemical A and B should come from different MAT files!')
    end

    % Call the functions that combine the MAT files for the two chemicals;
    % addition and independent action separately.
    par_out_ADD = plot_grid_add(filename);
    par_out_IND = plot_grid_ind(filename);
end

