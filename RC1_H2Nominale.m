clear all;
close all;
clc;
%% Data
M = 5;
m = 3;
L = 0.4;
J = m*L^2/3;
g = 9.81;
b = 0.5;
bmin = 0.7*b;
bmax = 1.3*b;
q = (M+m)*(J+m*L^2)-(m*L)^2;
p = q;
eps = 1*10^(-4);
%% Canonical State Space - Adapted (eps inserted)
Aeps1 = [0 1 0 0;
        0 0 1 0;
        0 0 0 1;
        bmin*g*m*L*eps/q, bmin*g*m*L/q, m*g*L*(M+m)/q,-bmin*(J+m*L^2)/q];
Aeps2 = [0 1 0 0;
        0 0 1 0;
        0 0 0 1;
        bmax*g*m*L*eps/q, bmax*g*m*L/q, m*g*L*(M+m)/q,-bmax*(J+m*L^2)/q];

Beps = [0;0;0;1];

Ceps = [-m*g*L/q, 0, (J+m*L^2)/q, 0;
        0, m*L*eps/q, m*L/q,0];

Deps = [0; 0];

Aeps = 0.5*(Aeps1+Aeps2);

%% H2 Optimal Control
% Select yhe matrices Q and R as follows
Q = Ceps'*Ceps;
R = eps^2;

P = care(Aeps,Beps,Q,R);
Ko = -inv(R)*(Beps'*P);
eigP = eig(P);
%% Plugging the controller to real system
An = [0 1 0 0;
         0 0 1 0;
         0 0 0 1;
         0 (b*g*m*L)/q m*g*L*(M+m)/q -b*(J+m*L^2)/q];
An1 = [0 1 0 0;
         0 0 1 0;
         0 0 0 1;
         0 (bmin*g*m*L)/q m*g*L*(M+m)/q -bmin*(J+m*L^2)/q];
An2 = [0 1 0 0;
     0 0 1 0;
     0 0 0 1;
     0 (bmax*g*m*L)/q m*g*L*(M+m)/q -bmax*(J+m*L^2)/q];
Bn = [0; 0; 0; 1];

Cn = [-m*g*L/q 0 (J+m*L^2)/q 0;
        0 0 m*L/q  0];
Dn = [0;0];

states = {'x1' 'x2' 'x3' 'x4'};
inputs = {'u'};
outputs = {'x [m]'; '\phi [rad]'};

ssControlled = ss(An+Bn*Ko,Bn, Cn, Dn,'statename',states,'inputname',...
    inputs,'outputname',outputs);
eig(An+Bn*Ko)

%% Plotting system response
 x0 = [-0.2544;0; 1.0705;0];
initial(ssControlled,x0)

%% H2 norm
H2Norm = trace(Bn'*P*Bn)

%% Quadratic Stability
% Stability Matrices
Aun1 = An1+Bn*Ko;
Aun2 = An2+Bn*Ko;

%Inequality
nstate = size(Aun1,1);
setlmis([]);    % Initialization of the LMI
P=lmivar(1, [nstate,1]);

% Subject function, LMI #1
% Aun1'P + PAun1 < 0
lmiterm([1 1 1 P], 1,Aun1 ,'s'); % LMI #1: Aun1'P + PAun1

% Subject function, LMI #2
% Aun2'P + P Aun2 < 0
lmiterm([2 1 1 P],1, Aun2,  's'); % LMI #1: Aun2'P + PAun2

% Subject function, LMI #3
% P>0
lmiterm([-3 1 1 P], 1, 1, 's'); % LMI #1: P>0

% Solving the LMI feasability problem
lmis = getlmis;
[tmin,xfeas] = feasp(lmis);
Popt = dec2mat(lmis,xfeas,P);

% Test for quadratic stability
eigPopt = eig(Popt)
quad1 = eig(Aun1'*Popt+Popt*Aun1)
quad2 = eig(Aun2'*Popt+Popt*Aun2)
