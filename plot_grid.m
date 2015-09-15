function plot_grid( grid )
% This is used to plot the grid that is being written for 2D cfd
% calculations

cornerlocs_x = grid.cornerlocs_x;
cornerlocs_y = grid.cornerlocs_y;

plot(cornerlocs_x, cornerlocs_y,'b.')

end

