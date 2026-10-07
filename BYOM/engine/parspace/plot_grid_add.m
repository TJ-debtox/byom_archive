function [par_out,error_flag] = plot_grid_add(filename)

% This is just to make the code backwards compatible (at least for a
% while). The function plot_grid_add has been renamed to
% plot_grid_add_guts, with the addition of plot_grid_add_debtox2019.

[par_out,error_flag] = plot_grid_add_guts(filename);