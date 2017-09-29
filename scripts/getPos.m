function out=getPos(~,input) %"~" sets default input if nothing is specified
global sunm G %global variables called from two_body_RK file

%% Initialise arrays
A2=zeros(6,1); %Final acceleration array 
A1=zeros(6,6); %zeros array to multiply with final acceleration
%Set values of acceleration in x,y and z to multiply with new acceleration
A1(1,4)=1; %corresponds to acceleration x
A1(2,5)=1; %corresponds to acceleration x
A1(3,6)=1; %corresponds to acceleration x

%% Calculate acceleration
x=input(1); %position x from finalPosition array
y=input(2); %position y from finalPosition array
z=input(3); %position z from finalPosition array
R=norm([x,y,z]); %Distance between planet and sun

unitvector=[x,y,z]/R; %Direction towards sun
A2(4)=-G*sunm/R^2*unitvector(1); %acceleration x for only planet
A2(5)=-G*sunm/R^2*unitvector(2); %acceleration y for only planet
A2(6)=-G*sunm/R^2*unitvector(3); %acceleration z for only planet

%% Output
%Total acceleration using the general "input" to calculate x,y and z
out=A1*input+A2; %multiply by A1 array to only 

end