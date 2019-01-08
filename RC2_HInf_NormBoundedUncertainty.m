close all;
clear all;
clc;

M = 5;
m = 3;
l = 0.4;
J = m*l^2/3;
g = 9.81;
b = 0.5;
q = (M+m)*(J+m*l^2)-(m*l)^2;

eps = 1*10^(-4);
nu = 1*10^(-10);

AhatEPS = [0 1 0 0;
        0 0 1 0;
        0 0 0 1;
        nu, nu, m*g*l*(M+m)/q, nu];
    
Delta_EPS = b*1.3;

L_EPS = [0;0;0;1];
N_EPS = [g*m*l*eps/q, g*m*l/q, 0, -1*(J+m*l^2)/q];

B_EPS = [0;0;0;1];

C_EPS = [-m*g*l/q, 0, (J+m*l^2)/q, 0;
    0, m*l*eps/q, m*l/q, 0];

D = [0; 0];

%% Real System
An = [0 1 0 0;
    0 0 1 0;
    0 0 0 1;
    0 (b*g*m*l)/q m*g*l*(M+m)/q -b*(J+m*l^2)/q];


AHatn = [0 1 0 0;
    0 0 1 0;
    0 0 0 1;
    0 nu m*g*l*(M+m)/q nu];
L_n = [0;0;0;1];
Delta_n = b;
N_n = [0 (g*m*l)/q 0 -1*(J+m*l^2)/q];

Bn = [0; 0; 0; 1];

Cn = [-m*g*l/q 0 (J+m*l^2)/q 0;
        0 0 m*l/q 0];
Dn = [0;0];
states = {'x1' 'x2' 'x3' 'x4'};
inputs = {'u'};
outputs = {'x [m]'; '\phi [rad]'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---------------------- Hinf Controller ------------------------------ %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AH2 = AhatEPS;
B2 = B_EPS;
B1 = L_EPS;
C1 = N_EPS;
C2 = C_EPS;
D12 = eps;

alpha = Delta_EPS;
gamma = inv(alpha);


Matr_B = [B2, B1];
Matr_R = [D12'*D12, 0; 0, inv(gamma^2-D12'*D12)];
Matr_S = [C1'*D12, C1'*D12];

P_HInf = care(AH2, Matr_B, C1'*C1, Matr_R, Matr_S);
K_o = -inv(D12'*D12)*(B2'*P_HInf+D12'*C1);

P_HInfv2 = care(AH2, Matr_B, C2'*C2, Matr_R, Matr_S);
K_ov2 = -inv(D12'*D12)*(B2'*P_HInfv2+D12'*C1);

%% Hinf Controller - Plugging the controller to real system

ssHInfv1 = ss(An+Bn*K_o,Bn, Cn, Dn,'statename',states,'inputname',...
    inputs,'outputname',outputs);
eig(An+Bn*K_o)

ssHInfv2 = ss(An+Bn*K_ov2,Bn, Cn, Dn,'statename',states,'inputname',...
    inputs,'outputname',outputs);
eig(An+Bn*K_ov2)

%% Hinf Controller - Plotting H2 system response
x0 = [-0.2544;0; 1.0705;0];
figure(1);
initial(ssHInfv1,x0);
figure(2);
initial(ssHInfv2,x0);

%% H2 norm
H2Tzw = trace(L_EPS'*P_HInf*L_EPS)
H2Norm = trace(Bn'*P_HInf*Bn)
H2Normv2 = trace(Bn'*P_HInfv2*Bn)

%% Quadratic Stability
s = tf('s');

quad = norm((N_n*inv(s*eye(4)-(An+Bn*K_o))*L_n),Inf)<inv(b)
quadv2 = norm((N_n*inv(s*eye(4)-(An+Bn*K_ov2))*L_n),Inf)<inv(b)

