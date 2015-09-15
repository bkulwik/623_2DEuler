function MakeSquareGridMatlab( initial_condition, filename )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

grid = struct;

%% Test Section

% Make sure test section is square
xmin = 0;
xmax = 8;
ymin = 0;
ymax = 8;

numcells_x = 8;
numcells_y = 8;


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
    
nodenumbers = [0:1:(length(x_nodes)*length(y_nodes))];

cellnum = 0;
cellnumbers = [];
cornerlocs_x = [];
cornerlocs_y = [];
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
        
        % Find adjacent cells
        if iX ~= 1
            cell_l = cellnum - 1;
        else
            cell_l = -1;
        end
        
        if iY ~= 1
            cell_b = cellnumbers(end) - (numcells_x+2);
        else
            cell_b = -1;
        end
        
        
        if iX ~= numcells_x + 2
            cell_r = cellnum + 1;
        else
            cell_r = -1;
        end
        
        if iY ~= numcells_y + 2
            cell_t = cellnum + (numcells_x + 2);
        else
            cell_t = -1;
        end
        
        adjacent_cells = [adjacent_cells; [cell_l cell_b cell_r cell_t]];
        
        if ((iX == 1) || (iX == length(x_cellcenter)) || (iY == 1) || (iY == length(y_cellcenter)))
            edge_cell = [edge_cell; 1];
        else
            edge_cell = [edge_cell; 0];
        end
        
        initial_conditions = [initial_conditions; initial_condition];
                
        cellnum = cellnum + 1;
    end
end




%% Write this all to a .bkcfd file
grid.cellnumber = cellnumbers;
grid.cornerlocs_x = cornerlocs_x;
grid.cornerlocs_y = cornerlocs_y;
grid.adjacent_cells = adjacent_cells;
grid.edge_cell = edge_cell;
grid.initial_conditions = initial_conditions;


struct2csv(grid,[filename '.bkcfd'])

end