close all;
clear all;
clc;

M = 5;
m = 3;
L = 0.4;
J = m*L^2/3;
g = 9.81;
b = 0.5;
bmin = b*0.7;
bmax = b*1.3;
q = (M+m)*(J+m*L^2)-(m*L)^2;

eps = 1*10^(-4);

A1_Eps = [0 1 0 0;
    0 0 1 0;
    0 0 0 1;
    bmin*g*m*L*eps/q, bmin*g*m*L/q, m*g*L*(M+m)/q,-bmin*(J+m*L^2)/q];
A2_Eps = [0 1 0 0;
    0 0 1 0;
    0 0 0 1;
    bmax*g*m*L*eps/q, bmax*g*m*L/q, m*g*L*(M+m)/q,-bmax*(J+m*L^2)/q];

B_Eps = [0;0;0;1];

C_Eps = [-m*g*L/q, 0, (J+m*L^2)/q, 0;
    0, m*L*eps/q, m*L/q,0];

D_Eps = [0; 0];

%% Real System
An = [0 1 0 0;
    0 0 1 0;
    0 0 0 1;
    0 (b*g*m*L)/q m*g*L*(M+m)/q -b*(J+m*L^2)/q];

A1n = [0 1 0 0;
    0 0 1 0;
    0 0 0 1;
    0 (bmin*g*m*L)/q m*g*L*(M+m)/q -bmin*(J+m*L^2)/q];
A2n = [0 1 0 0;
    0 0 1 0;
    0 0 0 1;
    0 (bmax*g*m*L)/q m*g*L*(M+m)/q -bmax*(J+m*L^2)/q];

Bn = [0; 0; 0; 1];

Cn = [-m*g*L/q 0 (J+m*L^2)/q 0;
        0 0 m*L/q 0];
Dn = [0;0];

states = {'x1' 'x2' 'x3' 'x4'};
inputs = {'u'};
outputs = {'x [m]'; '\phi [rad]'};

%% Control Matrices
D12 = eps;
Q = C_Eps'*C_Eps;
R = D12'*D12;

%% Enforcing simmetry
Q = (Q+Q.')/2;
Q_half = Q^(1/2);
Q_half=(Q_half+Q_half.')/2;

R = (R+R.')/2;
R_half = R^(1/2);
R_half=(R_half+R_half.')/2;
 
%% LMI configuration
% The system state space is:
% x_dot = Ai x + B u
% y = C x
%
% To solve the equation via LMI:
% [Z [C; 0]*S + [0; R1/2]W]>0
% A1*S+B*W+W'*B'+S*A1'+X0*X0'<0
% A2*S+B*W+W'*B'+S*A2'+X0*X0'<0
%
% The control is: u = W*S^(-1)*x
nstate = size(A1_Eps,1);

% Initialization of the LMI
setlmis([]);    
S=lmivar(1, [nstate,1]);
Z=lmivar(2, [1,1]);
W=lmivar(2, [1,nstate]);

% Subject function, LMI #1
% A1S + BW + SA1' + W'B' + BB' < 0
lmiterm([1 1 1 S], A1_Eps, 1, 's'); % LMI #1: A1S + SA1'
lmiterm([1 1 1 W], B_Eps, 1, 's'); % LMI #1: BW + W'B'
lmiterm([1 1 1 0],B_Eps*B_Eps'); % LMI #1: BB'

% Subject function, LMI #2
% A1S + BW + SA1' + W'B' + BB' < 0
lmiterm([2 1 1 S], A2_Eps, 1, 's'); % LMI #1: A2S + SA2'
lmiterm([2 1 1 W], B_Eps, 1, 's'); % LMI #1: BW + W'B'
lmiterm([2 1 1 0],B_Eps*B_Eps'); % LMI #1: BB'


% Subject function, LMI #3:
% [ Z Sqrt(R)W ]
% [ ] > 0
% [W'Sqrt(R) S ]
lmiterm([-3 1 1 Z], 1, 1); % LMI #2: Z
lmiterm([-3 2 1 -W], 1, R_half); % LMI #2: W'*sqrt(R)
lmiterm([-3 2 2 S], 1, 1); % LMI #2: S

% Subject function, LMI #4:
% S>0
lmiterm([-4 1 1 S], 1, 1, 's'); % LMI #4: S>0

% Create the LMI system
lmisys = getlmis;
n = decnbr(lmisys);
c = zeros(n,1);
for i=1:n
[Sj, Zj, Wj] = defcx(lmisys, i, S, Z, W);
c(i) =  trace(Zj)+trace(Q*Sj);
end

% Solving LMIs
options = [1e-20, 100, -1, 5, 1];
[copt, xopt] = mincx (lmisys, c, options);
% Results
Zopt = dec2mat(lmisys, xopt, Z);
Wopt = dec2mat(lmisys, xopt, W);
Sopt = dec2mat(lmisys, xopt, S);

% K stabilizing obtained via LMI
Klmi = Wopt*inv(Sopt);
eigSopt = eig(Sopt)

%% New System
Anew = An + Bn*Klmi;
eig(Anew)
ssNew = ss(Anew, Bn, Cn, Dn,'statename',states,'inputname',...
    inputs,'outputname',outputs);
figure(1)
 x0 = [-0.2544;0; 1.0705;0];
initial(ssNew,x0)

%% H2 norm
H2Norm = trace((Cn*Sopt+D12*Wopt)*inv(S)*(Cn*Sopt+D12*Wopt)') 

%% Quadratic Stability
Aun1 = A1n+Bn*Klmi;
Aun2 = A2n+Bn*Klmi;

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
Popt = dec2mat(lmis,xfeas,P)

% Test for quadratic stability
eigP = eig(Popt)
quad1 = eig(Aun1'*Popt+Popt*Aun1)
quad2 = eig(Aun2'*Popt+Popt*Aun2)