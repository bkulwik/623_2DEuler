function Plot_Solution( filename )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

close all
gamma = 1.4;

M = csvread(filename);
x_node_loc = M(:,1:4);
y_node_loc = M(:,5:8);
rho = M(:,9);
rhou = M(:,10);
rhov = M(:,11);
E = M(:,12);


for k = 1:size(x_node_loc,1)
    x_cellcenter(k) = mean(x_node_loc(k,:));
    y_cellcenter(k) = mean(y_node_loc(k,:));
end

x_cellcenter = unique(x_cellcenter);
y_cellcenter = unique(y_cellcenter);


[Xcell,Ycell] = meshgrid(x_cellcenter,y_cellcenter);

figure(1);
rho = convert_vec2mat(rho, size(x_cellcenter,2));
[~,h] = contourf(Xcell,Ycell,rho,200);
colormap(jet)
set(h,'EdgeColor','none');
colorbar;
xlabel('X position')
ylabel('Y position')
title('Density')

figure (2);
rhou = convert_vec2mat(rhou, size(x_cellcenter,2));
rhov = convert_vec2mat(rhov, size(x_cellcenter,2));
vel_mag = sqrt((rhou./rho).^2+(rhov./rho).^2);

[~,h] = contourf(Xcell,Ycell,vel_mag,200);
set(h,'EdgeColor','none');
colormap(jet)
colorbar;
xlabel('X position')
ylabel('Y position')
title('Velocity Magnitude')

figure(3);
P = (gamma-1).*(convert_vec2mat(E, size(x_cellcenter,2))-rho.*vel_mag.^2./2);
[~,h] = contourf(Xcell,Ycell,P,200);
set(h,'EdgeColor','none');
colormap(jet)
colorbar;
xlabel('X position')
ylabel('Y position')
title('Pressure')

end

function mat = convert_vec2mat(vector, len) 

rownum = 1;
while 1
    try
        mat(rownum,1:len) = vector(len*(rownum-1)+1:len*(rownum));
    catch
        break;
    end
    rownum = rownum + 1;
end
end