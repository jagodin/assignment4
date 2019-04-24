% netlist.m: 
%
% ELEC4700, Assignment 4
% Author: Jacob Godin
% Date: 2019/03/19
%--------------------------------------------------------------------------

%clear all;

global G C b; % define global variables

num_nodes = 6; % number of nodes in circuit

G = zeros(num_nodes,num_nodes); % Define G, 4 node circuit (do not include additional variables)
C = zeros(num_nodes,num_nodes); % Define C, 4 node circuit (do not include additional variables)
b = zeros(num_nodes,1); % Define b, 4 node circuit (do not include additional variables)

%--------------------------------------------------------------------------
% List of the components and values (netlist):
%--------------------------------------------------------------------------

R1 = 1;
R2 = 2;
R3 = 10;
R4 = 0.1;
R0 = 1000;
Cap = 0.25;
L = 0.2;
alpha = 100;

res(1,2,R1);
res(5,6,R4);
res(6,0,R0);
res(2,0,R2);
res(3,4,R3)
cap(1,2,Cap);

% Additional rows/columns stamps
vol(1,0,1); % Will change b vector on each iteration to input voltage vector
ind(2,3,L);
ccvs(5,0,4,0,alpha);
