close all;
clear all;
clc;

M = 5;
m = 3;
L = 0.4;
J = m*L^2/3;
g = 9.81;
b = 0.5;
q = (M+m)*(J+m*L^2)-(m*L)^2;

eps = 1*10^(-4);
nu = 1*10^(-10);

Ahateps = [0 1 0 0;
        0 0 1 0;
        0 0 0 1;
        nu, nu, m*g*L*(M+m)/q, nu];
    
Delta_eps = b*1.3;

L_eps = [0;0;0;1];
N_eps = [g*m*L*eps/q, g*m*L/q, 0, -1*(J+m*L^2)/q];

B_eps = [0;0;0;1];

C_eps = [-m*g*L/q, 0, (J+m*L^2)/q, 0;
    0, m*L*eps/q, m*L/q, 0];

D = [0; 0];

%% Real System
An = [  0 1 0 0;
        0 0 1 0;
        0 0 0 1;
        0, (b*g*m*L)/q, m*g*L*(M+m)/q, -b*(J+m*L^2)/q];


AHatn = [0 1 0 0;
        0 0 1 0;
        0 0 0 1;
        0, nu, m*g*L*(M+m)/q, nu];

L_n = [0;0;0;1];
Delta_n = b;
N_n = [0, (g*m*L)/q, 0, -1*(J+m*L^2)/q];

Bn = [0; 0; 0; 1];

Cn = [-m*g*L/q 0 (J+m*L^2)/q 0;
        0 0 m*L/q 0];
Dn = [0;0];
states = {'x1' 'x2' 'x3' 'x4'};
inputs = {'u'};
outputs = {'x [m]'; '\phi [rad]'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------------------- H2 Controller ----------------------------- %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AH2 = Ahateps;
B2 = B_eps;
B1 = L_eps;
C1 = N_eps;
C2 = C_eps;
D12 = eps;

P_H2 = care(AH2, B2, C1'*C1, D12'*D12, C1'*D12);
K_o = -inv(D12'*D12)*(B2'*P_H2+D12'*C1);

P_H2_v2 = care(AH2, B2, C2'*C2, D12'*D12, C1'*D12);
K_o_v2 = -inv(D12'*D12)*(B2'*P_H2_v2+D12'*C1);
%% H2 Controller - Plugging the controller to real system

ssH2v1 = ss(An+Bn*K_o,Bn, Cn, Dn,'statename',states,'inputname',...
    inputs,'outputname',outputs);
eig(An+Bn*K_o)

ssH2_v2 = ss(An+Bn*K_o_v2,Bn, Cn, Dn,'statename',states,'inputname',...
    inputs,'outputname',outputs);
eig(An+Bn*K_o_v2)
%% H2 Controller - Plotting H2 system response
x0 = [-0.2544;0; 1.0705;0];
figure(1);
initial(ssH2v1,x0)

figure(2);
initial(ssH2_v2,x0);

%% H2 Controller - H2 Norm
H2Tzw = trace(L_eps'*P_H2*L_eps)
H2Norm = trace(Bn'*P_H2*Bn)
H2Normv2 = trace(Bn'*P_H2_v2*Bn)

%% Quadratic Stability
s = tf('s');

quad = norm((N_n*inv(s*eye(4)-(An+Bn*K_o))*L_n),Inf)<inv(b)
quadv2 = norm((N_n*inv(s*eye(4)-(An+Bn*K_o_v2))*L_n),Inf)<inv(b)

% LMI double check for H2v2 controller - 
% must be (A+LDN)'P+P(A+LDN)<0 with P>0

% Stability Matrix
Aun = AHatn+Bn*K_o_v2+L_n*Delta_n*N_n;

%Inequality
nstate = size(Aun,1);
setlmis([]);    % Initialization of the LMI
P=lmivar(1, [nstate,1]);

% Subject function, LMI #1
% Aun1'P + PAun1 < 0
lmiterm([1 1 1 P], 1,Aun ,'s'); % LMI #1: Aun1'P + PAun1


% Subject function, LMI #3
% P>0
lmiterm([-2 1 1 P], 1, 1, 's'); % LMI #1: P>0

% Solving the LMI feasability problem
lmis = getlmis;
[tmin,xfeas] = feasp(lmis);
Popt = dec2mat(lmis,xfeas,P)

% Test for quadratic stability
eigPopt = eig(Popt)
quadv2INEQ = eig(Aun'*Popt+Popt*Aun)