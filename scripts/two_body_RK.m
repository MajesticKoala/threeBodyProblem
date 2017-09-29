clear all
close all
clc
hold on

dt=1;
finaltime=1000; %Total time in days

%% Graphing
range = 0.02;
axis([-range range -range range -1e-6 1e-6]);
title('Orbits');
xlabel('displacement x AU');
ylabel('displacement y AU');
whitebg('black');
view(-20, 20);
grid on

%% Variables
global sunm G %Global constants so they can be used in function
sunm = 1.9891e30; %kg
G = 1.4879e-34;   %kg m AU

%% Initial conditions
%sun positioned at origin
x=0;
y=0;
z=0;
vx=0;
vy=0;
vz=0;
SunInit=[x,y,z,vx,vy,vz]; %Build array of planet components

%Earth
x=-0.78297;     %x position
y=0.6009;       %y position
z=-1.3779E-05;  %z position
vx=-1.075E-02;  %x velocity
vy=-1.3705E-02; %y velocity
vz=3.974E-08;   %z velocity
EarthInit=[x,y,z,vx,vy,vz]; %Build array of planet components

%% Call function
%Call ode45 function on getPos() which calculates the new position
[attime, finalposition]=ode45(@getPos,(0:dt:finaltime),(EarthInit), (SunInit));

% Plot
[X,Y] = meshgrid(-range:0.005:range); %Set meshgrid to display xy plane
Z=X*0;                                %Z value = 0 so flat plane
mesh(X,Y,Z);                          %Displays Meshgrid
alpha(0)                              %Transparent
plot(0, 0, 'white o','MarkerFaceColor','white', 'MarkerSize', 10) %plots the sun at centre of system

%%Animating
for i = 1:finaltime
   %Animate every point in finalposition array
   s = plot3(finalposition(i,4),finalposition(i,5), finalposition(i,6), 'red o','MarkerFaceColor','red', 'MarkerSize', 5); 
   pause(0.01);
   delete(s);
end

%Plots the total planet orbit
plot3(finalposition(:,4), finalposition(:,5), finalposition(:,6), 'g');
