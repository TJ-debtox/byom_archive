function [par_out,error_flag] = plot_grid_ind(filename)

% This is just to make the code backwards compatible (at least for a
% while). The function plot_grid_ind has been renamed to
% plot_grid_ind_guts, with the addition of plot_grid_ind_debtox2019.

[par_out,error_flag] = plot_grid_ind_guts(filename);