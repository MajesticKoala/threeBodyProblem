clear all
close all
clc
%% Initialise variables
m1 = 1;       % mass 1
m2 = 0.1;     % mass 2
M = m1 + m2;  % total mass


T = 2*pi * sqrt(1 / M);             % period of orbit
omega = 2*pi/T;                     % angular velocity based on orbit

[X,Y] = meshgrid(-1.5:0.01:1.5);    % meshgrid sets out entire plane
U = calcPotential(m1, m2, X, Y);    % total potential energy to be plotted

%% Plotting
surf(X,Y,U, 'EdgeColor', 'none');   % plots total potential energy for every point on plane
colormap('hot');
caxis([-4 0]);                      % sets the colour plot to a small range so Lagrange points can be seen
title("Lagrange Point map");

axis([-1 1.5 -1 1]);                % axis limits
zlim([-4 0]);                       % essentially sets a cutoff in z (-4, -2)
view(90,90)                         % top down view

%% Orbiting camera
for i = 1:1000                      %Rotates camera around the plot
    surf(X,Y,U, 'EdgeColor', 'none');
    axis([-1 1 -1 1]);
    caxis([-4 -2]);
    zlim([-4 -2]); 
   view(90-i, 50);
   pause(0);
end

function U = calcPotential(m1, m2, x, y)

    K = 1.7.*sqrt(m1+m2);           % angular velocity

    x1 = -m2/(m1+m2);               % x position of m1, y = 0
    x2 = 1 + x1;                    % x position of m2, y = 0

    r1 = sqrt((x-x1).^2+y.^2);      % distance from m1 to point
    r2 = sqrt((x-x2).^2+y.^2);      % distance from m2 to point

    U = -K^2/2*(x.^2+y.^2) - m1./r1 - m2./r2; % total potential at that point  
end
