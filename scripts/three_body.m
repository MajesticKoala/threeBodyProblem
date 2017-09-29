clear all
clc
close all
hold on

% Graphing
range = 2;
axis([-range range -range range]); %Sets a fixed axis
title('Orbits');
xlabel('displacement x AU');
ylabel('displacement y AU');
whitebg('black');

txtTime = text(-range+0.2, range-0.4, 's', 'fontsize', 16); %Displays a timer in top left of screen

%% Constants
global dt G;            %Set global variables so they can be used inside functions
G = 1.4879e-34;         % AU m kg
mass_sun = 1.9891e30;   %kg
mass_earth = 5.972e24; 
mass_moon = 7.348e22;

dt = 3;                 %time step
finaltime = 10000;      %final time
t = (1:dt:finaltime);   %Initialises time array


%% Initialise velocity, acceleration and distance arrays
%1st body
vx1 = t*0;  %x velocity (same size as time array)
vy1= t*0;   %y velocity
ax1= t*0;   %x acceleration
ay1= t*0;   %y acceleration
dx1= t*0;   %x position
dy1= t*0;   %y position

%2nd body
vx2= t*0;
vy2= t*0;
ax2= t*0;
ay2= t*0;
dx2= t*0;
dy2= t*0;

%3rd body
vx3= t*0;
vy3= t*0;
ax3= t*0;
ay3= t*0;
dx3= t*0;
dy3= t*0;

%% Initial Conditions
%%Initial Position
                    %Circle Data                  %Figure 8 Data
dx1(1)=-1.7;             %-1.7;                        %-0.97000436;
dy1(1)=0;             %0;                           %0.24308753;

dx2(1)=1.7;             %1.7;                         %0.97000436;
dy2(1)=0;             %0;                           %-0.24308753;

dx3(1)=0;
dy3(1)=0;

%Initial Velocities (positive is left, down)
%(earth 0.01732) AU / day
%Constants used for figure 8 orbits %AU/day 
const1 = 0.0008794;
const2 = 0.0008117;

                    %Circle Data                %Figure 8 data
vx1(1)=0;             %0;                         %-const1/2;
vy1(1)=-0.0008;             %-0.0008;                   %-const2/2;

vx2(1)=0;             %0;                         %-const1/2;
vy2(1)=0.0008;             %0.0008;                    %-const2/2; 

vx3(1)=0;             %0;                         %const1;
vy3(1)=0;             %0;                         %const2;

%initial masses 
m1=mass_earth.*1000; %mass for mass1
m2=mass_earth.*1000; %mass for mass2
m3=mass_earth.*1000; %mass for mass3

%Energy
kinetic1 = t*0; %Kinetic energy for mass1
kinetic2 = t*0; %Kinetic energy for mass2
kinetic3 = t*0; %Kinetic energy for mass3
kinTotal = t*0; %Total Kinetic energy

potential1 = t*0; %Potential energy for mass1
potential2 = t*0; %Potential energy for mass2
potential3 = t*0; %Potential energy for mass3
potTotal = t*0;   %Total Potential energy

energyTotal = t*0;%Total energy of system

kinetic1(1) = 0.5.*m1.*(vx1(1).^2 +vy1(1).^2); %Initialise Kinetic energy for mass1
kinetic2(1) = 0.5.*m2.*(vx2(1).^2 +vy2(1).^2); %Initialise Kinetic energy for mass2
kinetic3(1) = 0.5.*m3.*(vx3(1).^2 +vy3(1).^2); %Initialise Kinetic energy for mass3

potential1(1) = (G*(m1*m2)./(getDistance(dx1(1), dy1(1), dx2(1), dy2(1)))) + (G*(m1*m3)./(getDistance(dx1(1), dy1(1), dx3(1), dy3(1)))); %Initialise Potential energy for mass1
potential2(1) = (G*(m2*m1)./(getDistance(dx1(1), dy1(1), dx2(1), dy2(1)))) + (G*(m2*m3)./(getDistance(dx2(1), dy2(1), dx3(1), dy3(1)))); %Initialise Potential energy for mass2
potential3(1) = (G*(m3*m1)./(getDistance(dx1(1), dy1(1), dx3(1), dy3(1)))) + (G*(m3*m2)./(getDistance(dx2(1), dy2(1), dx3(1), dy3(1)))); %Initialise Potential energy for mass3


%% Main Loop
tic
for i = 1:length(t)-1
    
%Energy%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        %calculate Kinetic energy using velocity and mass
        kinetic1(i) = 0.5.*m1.*(vx1(i).^2 +vy1(i).^2); 
        kinetic2(i) = 0.5.*m2.*(vx2(i).^2 +vy2(i).^2);
        kinetic3(i) = 0.5.*m3.*(vx3(i).^2 +vy3(i).^2);
        kinTotal(i) = kinetic1(i)+kinetic2(i)+kinetic3(i); %Update total kinetic energy as sum of all components

        %calculate Potential energy for each mass with respect to the other two
        potential1(i) = (G*(m1*m2)./(getDistance(dx1(i), dy1(i), dx2(i), dy2(i)))) + (G*(m1*m3)./(getDistance(dx1(i), dy1(i), dx3(i), dy3(i)))); 
        potential2(i) = (G*(m2*m1)./(getDistance(dx1(i), dy1(i), dx2(i), dy2(i)))) + (G*(m2*m3)./(getDistance(dx2(i), dy2(i), dx3(i), dy3(i))));
        potential3(i) = (G*(m3*m1)./(getDistance(dx1(i), dy1(i), dx3(i), dy3(i)))) + (G*(m3*m2)./(getDistance(dx2(i), dy2(i), dx3(i), dy3(i)))); 
        potTotal(i) = potential1(i) + potential2(i)+ potential3(i); %Update total potential energy as sum of all components

        energyTotal(i) = potTotal(i) + kinTotal(i); %Sum of potential energy and kinetic energy
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %mass1 (Largest if different masses)
    %Calculate force and angle with respect to other two masses using
    %functions getForce() and getAngle()
    force12 = getForce(dx1(i), dy1(i), dx2(i), dy2(i), m1, m2);
    force13 = getForce(dx1(i), dy1(i), dx3(i), dy3(i), m3, m1);
    angle12 = getAngle(dx1(i), dy1(i), dx2(i), dy2(i));
    angle13 = getAngle(dx1(i), dy1(i), dx3(i), dy3(i));
    
    %Acceleration using force and angle using function getAcc()
    ax1(i) = getAccx(force12, angle12, m1)+getAccx(force13, angle13, m1);
    ay1(i) = getAccy(force12, angle12, m1)+getAccy(force13, angle13, m1);
    %Determine velocity from acceleration using function getVel()
    vx1(i+1) = getVelx(ax1(i), vx1(i));
    vy1(i+1) = getVely(ay1(i), vy1(i));
    %Calculate new position using function getDis()
    dx1(i+1) = getDisx(dx1(i), vx1(i), ax1(i));
    dy1(i+1) = getDisy(dy1(i), vy1(i), ay1(i));
    
    %plot the new position
    p1 = plot(dx1(i), dy1(i),'red o','MarkerFaceColor','red', 'MarkerSize', 5);
    
    %mass 2 (Second largest if different masses)
    force21 = getForce(dx1(i), dy1(i), dx2(i), dy2(i), m1, m2);
    force23 = getForce(dx2(i), dy2(i), dx3(i), dy3(i), m3, m2);
    angle21 = getAngle(dx2(i), dy2(i), dx1(i), dy1(i));
    angle23 = getAngle(dx2(i), dy2(i), dx3(i), dy3(i));
    
    ax2(i) = getAccx(force21, angle21, m2)+getAccx(force23, angle23, m2);
    ay2(i) = getAccy(force21, angle21, m2)+getAccy(force23, angle23, m2);
    vx2(i+1) = getVelx(ax2(i), vx2(i));
    vy2(i+1) = getVely(ay2(i), vy2(i));
    dx2(i+1) = getDisx(dx2(i), vx2(i), ax2(i));
    dy2(i+1) = getDisy(dy2(i), vy2(i), ay2(i));
    p2 = plot(dx2(i), dy2(i),'blue o','MarkerFaceColor','blue', 'MarkerSize', 3);

    %mass 3 (Smallest if different masses)
    force31 = getForce(dx1(i), dy1(i), dx3(i), dy3(i), m1, m3);
    force32 = getForce(dx2(i), dy2(i), dx3(i), dy3(i), m2, m3);
    angle31 = getAngle(dx3(i), dy3(i), dx1(i), dy1(i));
    angle32 = getAngle(dx3(i), dy3(i), dx2(i), dy2(i));
    
    ax3(i) = getAccx(force31, angle31, m3)+getAccx(force32, angle32, m3);
    ay3(i) = getAccy(force31, angle31, m3)+getAccy(force32, angle32, m3);
    vx3(i+1) = getVelx(ax3(i), vx3(i));
    vy3(i+1) = getVely(ay3(i), vy3(i));
    dx3(i+1) = getDisx(dx3(i), vx3(i), ax3(i));
    dy3(i+1) = getDisy(dy3(i), vy3(i), ay3(i));
    p3 = plot(dx3(i), dy3(i),'yellow o','MarkerFaceColor','yellow', 'MarkerSize', 2);
    
    pause(0); %Allows the graph to update even though 0s pause
    %Deletes all points to make an animation
    delete(p1);
    delete(p2);
    delete(p3);
    
    timeSimulation = i.*dt; %How many days the simulation has been running
    set(txtTime , 'string', [num2str(round(timeSimulation/365.25)) ' yrs']); %Displays how long the simulation has been running in years

end
toc %Simulation run time
plot(dx1, dy1, dx2, dy2, dx3, dy3); %Plots path of planets

%% Energy Graphs
%kinetic energy of each planet
figure                                          %New figure
whitebg('white');                               %White background
plot(t, kinetic1, t, kinetic2, t, kinetic3);    %Plots t vs kinetic energy
title('Kinetic energy of each planet'); 
xlabel('time');
ylabel('Kinetic Energy kg AU^2 days^-2');       %Units of Joule kg AU^2 days^-2
legend('mass1', 'mass2', 'mass3');

%total kinetic energy of all planets
figure
plot(t, kinTotal);
title('Total Kinetic Energy');
xlabel('time');
ylabel('Kinetic Energy kg AU^2 days^-2');

%potential energy of each planet
figure
plot(t, potential1, t, potential2, t, potential3); 
title('Potential energy of each planet');
xlabel('time');
ylabel('Potential Energy kg AU^2 days^-2');
legend("mass1", 'mass2', 'mass3');

%total potential energy of all planets
figure
plot(t, potTotal); 
title('Total Potential Energy');
xlabel('time');
ylabel('Potential Energy kg AU^2 days^-2');

%total energy of entire system
figure
plot(t, energyTotal); 
title('Total Energy');
xlabel('time');
ylabel('Total Energy kg AU^2 days^-2');

%Plots total kinetic, potential and the sum of both on one graph
figure
plot(t, potTotal, t, kinTotal, t, energyTotal);
title('Total Energy');
xlabel('time');
ylabel('Energy kg AU^2 days^-2');
legend("Potential", 'Kinetic', 'Total Energy');

%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get Force
function force = getForce(dx, dy, dx2, dy2, m, m2)
global G; %Calls the global variable
    force = (G.*m.*m2)./((getDistance(dx, dy, dx2, dy2).^2));
end

%get Distance
function distance = getDistance(dx, dy, dx2, dy2)
    distance = sqrt((dx-dx2).^2+(dy-dy2).^2);
end

%get Angle
function angle = getAngle(dx1, dy1, dx2, dy2)
    angle = atan2(dy1-dy2,dx1-dx2); %Atan2 return the angle between two points
end

%get Acceleration
function accx = getAccx(force, angle, mass)
    accx = (force.*cos(angle))./mass;
end
function accy = getAccy(force, angle, mass)
    accy = (force.*sin(angle))./mass;
end

%Get Velocity
function velx = getVelx(accxi, velxi)
global dt;
    velx = velxi + accxi.*dt;
end
function vely = getVely(accyi, velyi)
global dt;
    vely = velyi + accyi.*dt;
end

%get Distance
function disx = getDisx(disxi, velx, accx)
global dt;
    disx = disxi - velx.*dt + 0.5.*accx.*dt.^2;
end
function disy = getDisy(disyi, vely, accy)
global dt;
    disy = disyi - vely.*dt + 0.5.*accy.*dt.^2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%