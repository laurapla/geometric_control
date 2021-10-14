% University of California, Irvine - Spring 2021
% Laura Pla Olea - lplaolea@uci.edu

clc; close all; clear all;

%% State Space

syms x1 x2 x3 x4 x9 x10 x11 real;
Q = [x1 x2 x3 x4 x9 x10 x11].'; % Coordinates of the system

syms a da dh real; % Angle of attack, pitching velocity, plunging velocity
syms CNa beta A1 A2 b1 b2 M U c real; % Circulatory component constants
syms Ka Kq Ti real; % Noncirculatory component constants
syms Tp Tf Tv; % Nonlinear flow components
syms fp(x9); % Steady lift coefficient (is a function of alpha and plunging)
assume(fp(x9),'real');

% Pitch rate
q = da/(U*c);

% Drift vector
f = [-b1*beta^2*(2*U/c)*x1;
    -b2*beta^2*(2*U/c)*x2;
    -x3/(Ka*Ti);
    -x4/(Kq*Ti);
    CNa*beta^2*(2*U/c)*A1*b1/Tp*x1+CNa*beta^2*(2*U/c)*A2*b2/Tp*x2-4/(M*Ka*Ti*Tp)*x3-1/(M*Kq*Ti*Tp)*x4-x9/Tp;
    -x10/Tf+fp/Tf;
    -x11/Tv-CNa*beta^4*(2*U/c)^2*(A1*b1^2*x1+A2*b2^2*x2)*(1-((1+sqrt(x10))/2)^2)+CNa*beta^2*(2*U/c)*(A1*b1*x1+A2*b2*x2)*(-x10+fp)/(4*Tf)*(1/sqrt(x10)+1)];
f = formula(f);

ga = [1 1 1 0 4/(M*Tp) 0 CNa*beta^2*(2*U/c)*(A1*b1+A2*b2)*(1-((1+sqrt(x10))/2)^2)].'; % Control vector of the angle of attack
gq = [0.5 0.5 0 1 1/(M*Tp) 0 CNa*beta^2*(U/c)*(A1*b1+A2*b2)*(1-((1+sqrt(x10))/2)^2)].'; % Control vector of the plunging

%% Lie Brackets and Symmetric Products

% Symmetric products
gaga = symmetric_product(Q,f,ga,ga); % bad symmetric product (potential obstruction to STLC)
gagq = symmetric_product(Q,f,ga,gq); % good symmetric product
gqgq = symmetric_product(Q,f,gq,gq); % bad symmetric product (potential obstruction to STLC)

%% Lift coefficient

Cnia = -4/M*x3/(Ka*Ti)+4/M*(a+atan(dh/U));
Cniq = -1/M*x4/(Kq*Ti)+q/M;

Cnf = CNa*beta^2*(2*U/c)*(A1*b1*x1+A2*b2*x2)*((1+sqrt(x10))/2)^2;
Cnv = x11;
Cni = Cnia+Cniq;
H = Cnf+Cnv+Cni; % Lift (output of the system)

% Lie derivatives (rate of change of H along the dynamics given by gxgx)
LgagaH = Lie_Derivative(Q,gaga,H);
LgagqH = Lie_Derivative(Q,gagq,H);
LgqgqH = Lie_Derivative(Q,gqgq,H);

% % Fliess Functional expansion
% G = [f ga gh]; % Matrix of the drift and control vectors
% LH = Fliess_coefficients(Q,G,H,2);
% 
% % %% Drag
% % 
% % S = sqrt(2)*(2*(xc+GammaQS/2)/(pi*c)-c*da/2)/2;
% % D = H*tan(a)-pi*rho*c*S^2/2; % Drag (another possible output of the system)
% % 
% % % Lie derivatives (rate of change of Dmaybe along the dynamics given by gxgx)
% % LgagaDmaybe = Lie_Derivative(Q,gaga,D);
% % LgaghDmaybe = Lie_Derivative(Q,gagh,D);
% % LghghDmaybe = Lie_Derivative(Q,ghgh,D);
% % 
% % % Fliess functional expansion
% % LDmaybe = Fliess_coefficients(Q,G,D,2);