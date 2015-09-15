function MakeGridMatlab_definewholegridbutusepart( initial_condition_global, initial_condition_postshock, filename )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

grid = struct;

%% Test Section

% Make sure test section is square
xmin = 0;
xmax = 100;
ymin = -50;
ymax = 50;

numcells_x = 120;
numcells_y = 120;

% this is all in meters, not in cell-widths
tunnel_width = 10;
tunnel_length = 50;

sidetunnel_height = 15; %Must have more than two cells!
sidetunnel_lengthback = 20;
sidetunnel_lengthforward = 20;
sidetunnel_width = 5;
sidetunnel_xstart = 25;

grid_plotting = 1;


% Calculation for setup
dx = (xmax-xmin)/numcells_x;
dy = (ymax-ymin)/numcells_y;

x_nodes = linspace(xmin-dx,xmax+dx,numcells_x+3);
y_nodes = linspace(ymin-dy,ymax+dy,numcells_y+3);
[X_nodes, Y_nodes] = meshgrid(x_nodes,y_nodes);
X_nodes = X_nodes';
Y_nodes = Y_nodes';

x_cellcenter = linspace(xmin+dx/2, xmax-dx/2, numcells_x+2);
y_cellcenter = linspace(ymin+dy/2, ymax-dy/2, numcells_y+2);
    
cellnum = 0;
cellnumbers = [];
cornerlocs_x = [];
cornerlocs_y = [];
centerloc = [];
adjacent_cells = [];
initial_conditions = [];
edge_cell = [];


for iY = 1:length(y_cellcenter)
    for iX = 1:length(x_cellcenter)
                
        % Find cell number
        cellnumbers = [cellnumbers; cellnum];
        
        % Find node numbers on corners
        % Corner goes in order form BL, BR, TR, TL
        corner(1) = (iY-1)*(length(x_nodes)) + iX;
        corner(2) = (iY-1)*(length(x_nodes)) + iX + 1;
        corner(3) = (iY)*(length(x_nodes)) + iX + 1;
        corner(4) = (iY)*(length(x_nodes)) + iX;
        
        cornerlocs_x = [cornerlocs_x; X_nodes(corner)];
        cornerlocs_y = [cornerlocs_y; Y_nodes(corner)];
        centerloc = [centerloc; [mean(X_nodes(corner)), mean(Y_nodes(corner))]];
        
        % Find adjacent cells
        if iX ~= 1
            cell_l = cellnum - 1;
        else
            cell_l = -2;
        end
        
        if iY ~= 1
            cell_b = cellnumbers(end) - (numcells_x+2);
        else
            cell_b = -2;
        end
        
        
        if iX ~= numcells_x + 2
            cell_r = cellnum + 1;
        else
            cell_r = -2;
        end
        
        if iY ~= numcells_y + 2
            cell_t = cellnum + (numcells_x + 2);
        else
            cell_t = -2;
        end
        
        adjacent_cells = [adjacent_cells; [cell_l cell_b cell_r cell_t]];

        cellnum = cellnum + 1;
    end
end



%% Figure out which parts of the solution I actually care about

sidetunnel_xend = sidetunnel_xstart + sidetunnel_width;

tunnelcells = [];
expansioncells = [];
sidetunnelcells_top = [];
boundarycells = [];
adjacent_cells_index = adjacent_cells+1;

for it = 1:length(centerloc)
    if      ((centerloc(it,1) < tunnel_length) && ...
            (centerloc(it,1) > 0) && ...
            (centerloc(it,2) < tunnel_width/2) && ...plot_grid
            (centerloc(it,2) > -tunnel_width/2))
        tunnelcells = [tunnelcells; cellnumbers(it)];
        boundarycells = add_to_boundarycells(it,adjacent_cells, boundarycells);
    elseif  (centerloc(it,1) > tunnel_length) &&  (centerloc(it,1) > 0) &&...
            (centerloc(it,1) < xmax) &&...
            (centerloc(it,2) > -ymax) &&...
            (centerloc(it,2) < ymax)
        expansioncells = [expansioncells; cellnumbers(it)];   
        boundarycells = add_to_boundarycells(it,adjacent_cells, boundarycells);
    elseif  (centerloc(it,2) > tunnel_width/2) && ...
            (centerloc(it,2) < tunnel_width/2+sidetunnel_height) && ...
            (((centerloc(it,1) > sidetunnel_xstart) && ...
            (centerloc(it,1) < sidetunnel_xend)) || ...
            ((centerloc(it,1) > sidetunnel_xstart-sidetunnel_lengthback) && ...
            (centerloc(it,1) < sidetunnel_xstart + sidetunnel_lengthforward) && ...
            (centerloc(it,2) > tunnel_width/2+(sidetunnel_height-sidetunnel_width)) && ...
            (centerloc(it,2) < tunnel_width/2+(sidetunnel_height+sidetunnel_width))))
        sidetunnelcells_top = [sidetunnelcells_top; cellnumbers(it)];
        boundarycells = add_to_boundarycells(it,adjacent_cells, boundarycells);
    end
end
% in C++ notation
gridcells = sort([tunnelcells; expansioncells; sidetunnelcells_top]);
% in MATLAB notation
gridcells_index = gridcells+1;

%% Assign "edge_cell" to the newly-edged cells
boundarycells = unique(boundarycells);
boundarycells = boundarycells(~ismember(boundarycells,gridcells));

boundarycells_index = boundarycells + 1;

all_used_cells = sort([boundarycells; gridcells]);
all_used_cells_index = all_used_cells + 1;

% edge_state = 1: standard wall ghost cell
% edge_state = 2: farfield zero gradient ghost cell
% edge_state = 0: interior cell

edge_state = zeros(length(all_used_cells),1);

for it = 1:length(all_used_cells_index)
    if (centerloc(all_used_cells_index(it),1) > xmax) || ...
           (centerloc(all_used_cells_index(it),2) > ymax) || (centerloc(all_used_cells_index(it),2) < ymin)
        edge_state(it) = 2;
    elseif (ismember(all_used_cells_index(it), boundarycells_index))
        edge_state(it) = 1;       
    end
end


%% Define initial conditions
initialconditions = [initial_condition_global(1).*ones(length(all_used_cells),1), ...
    initial_condition_global(2).*ones(length(all_used_cells),1), ...
    initial_condition_global(3).*ones(length(all_used_cells),1), ...
    initial_condition_global(4).*ones(length(all_used_cells),1)];

initialconditions(find(centerloc(all_used_cells_index,1) < 0),1) = initial_condition_postshock(1);
initialconditions(find(centerloc(all_used_cells_index,1) < 0),2) = initial_condition_postshock(2);
initialconditions(find(centerloc(all_used_cells_index,1) < 0),3) = initial_condition_postshock(3);
initialconditions(find(centerloc(all_used_cells_index,1) < 0),4) = initial_condition_postshock(4);

initialconditions(find(centerloc(all_used_cells_index,1) < 0)+1,1) = initial_condition_postshock(1);
initialconditions(find(centerloc(all_used_cells_index,1) < 0)+1,2) = initial_condition_postshock(2);
initialconditions(find(centerloc(all_used_cells_index,1) < 0)+1,3) = initial_condition_postshock(3);
initialconditions(find(centerloc(all_used_cells_index,1) < 0)+1,4) = initial_condition_postshock(4);


% [all_used_cells, edge_state, initialconditions]

%% Write this all to a .bkcfd file

grid.cellnumber = all_used_cells;
grid.cornerlocs_x = cornerlocs_x(all_used_cells_index,:);
grid.cornerlocs_y = cornerlocs_y(all_used_cells_index,:);
grid.adjacent_cells = adjacent_cells(all_used_cells_index,:);
grid.edge_state = edge_state;
grid.initial_conditions = initialconditions;


struct2csv(grid,[filename '.bkcfd'])

%% Plot the grid to see a visual depiction of what was made

if grid_plotting
    plot_grid( grid )
end

end



%%
function boundarycells = add_to_boundarycells(index, adjacent_cells, boundarycells)

    boundarycells = [boundarycells; adjacent_cells(index,:) ];

end


