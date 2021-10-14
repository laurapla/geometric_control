% University of California, Irvine - Summer 2021
% Laura Pla Olea - lplaolea@uci.edu
tic
clc; close all; clear;

%% State Space

syms x1 x2 x3 x4 x9 x10 x11 real;
syms a da dh real; % Angle of attack, pitching velocity, plunging velocity
Q = [x1 x2 x3 x4 x9 x10 x11 a da dh].'; % Coordinates of the system


syms CNa beta A1 A2 b1 b2 M V c b real; % Circulatory component constants
syms Ka Kq Ti real; % Noncirculatory component constants
syms Tp Tf Tv; % Nonlinear flow components
syms fp(x1,x2,x3,x4,x9,x10,a,da,dh); % Steady lift coefficient (is a function of alpha and plunging)
assume(fp(x1,x2,x3,x4,x9,x10,a,da,dh),'real');
syms dCv(x1,x2,x3,x4,x9,x10,a,da,dh);
assume(dCv(x1,x2,x3,x4,x9,x10,a,da,dh),'real');

% Drift vector
f = [-b1*beta^2*(2*V/c)*x1+a+atan(dh/V)+da*c/(2*V);
    -b2*beta^2*(2*V/c)*x2+a+atan(dh/V)+da*c/(2*V);
    -x3/(Ka*Ti)+a+atan(dh/V);
    -x4/(Kq*Ti)+da*c/V;
    (CNa*beta^2*(2*V/c)*(A1*b1*x1+A2*b2*x2)-4/(M*Ka*Ti)*x3-1/(M*Kq*Ti)*x4-x9+4/M*(a+atan(dh/V))+da*c/V/M)*2*V/c/Tp;
    (-x10/Tf+fp/Tf)*2*V/c;
    (-x11/Tv+dCv/Tv)*2*V/c;
    da;
    0;
    0];
f = formula(f);

ga = [0 0 0 0 0 0 0 0 1 0].'; % Control vector of the pitching motion
gh = [0 0 0 0 0 0 0 0 0 1].'; % Control vector of the plunging motion

% %% Lie Brackets and Symmetric Products
% 
% % Symmetric products
% gaga = symmetric_product(Q,f,ga,ga); % bad symmetric product (potential obstruction to STLC)
% gagh = symmetric_product(Q,f,ga,gh); % good symmetric product
% ghgh = symmetric_product(Q,f,gh,gh); % bad symmetric product (potential obstruction to STLC)
% 
%% Lift and drag coefficients

q = da*c/V; % Pitch rate

% Normal force coefficients
Cnia = -4/M*x3/(Ka*Ti)+4/M*(a+atan(dh/V)); % Non-circulatory load due to the angle of attack
Cniq = -1/M*x4/(Kq*Ti)+q/M; % Non-circulatory load due to the pitch rate
Cni = Cnia+Cniq; % Total non-circulatory (impulsive) normal force coefficient

Cnf = CNa*beta^2*(2*V/c)*(A1*b1*x1+A2*b2*x2)*((1+sqrt(x10))/2)^2; % Circulatory load due to TE separation
Cnv = x11; % Circulatory load due to LE separation (vortex)

Cn = Cnf+Cnv+Cni; % Normal force coefficient

% Chord force coefficient
aE = (2*V/c)*beta^2*(A1*b1*x1+A2*b2*x2); % Effective angle of attack
syms eta real; % Parameter that accounts for the viscous effects
Ccp = CNa*aE^2; % Chord force under attached flow conditions
Ccf = eta*CNa*aE^2*sqrt(x10); % Chord force coefficient under TE separation
Cc = Ccf; % Chord force coefficient

% Total lift and drag coefficients
Cl = Cn*cos(a)+Cc*sin(a); % Lift coefficient
Cd = Cn*sin(a)-Cc*cos(a); % Drag coefficient
% 
% H = Cl;
% D = Cd;
% 
% % Lie derivatives (rate of change of H along the dynamics given by gxgx)
% LgagaH = Lie_Derivative(Q,gaga,H);
% LgaghH = Lie_Derivative(Q,gagh,H);
% LghghH = Lie_Derivative(Q,ghgh,H);
% 
% % Lie derivatives (rate of change of D along the dynamics given by gxgx)
% LgagaD = Lie_Derivative(Q,gaga,D);
% LgaghD = Lie_Derivative(Q,gagh,D);
% LghghD = Lie_Derivative(Q,ghgh,D);
% 
% %% Pullback
% 
% n = 2; % Pullback order
% 
syms w A_alpha H phi real; % Constants of the inputs
% syms t real; % Independent variable
% 
% G = [ga gh];
% 
% % Inputs
% u_alpha = w^2*A_alpha*cos(w*t);
% u_h = w^2*H*b*cos(w*t+phi);
% Us = [u_alpha u_h];
% 
% % Pullback of f along the flow of G
% F = simplify(pullback(Q,f,G,Us,t,n),'IgnoreAnalyticConstraints',true);
% Gvector_averaged = simplify(w/(2*pi)*int(F-f,t,0,2*pi/w));
% 

%% Series expansion (averaged values)

% Amplitude and phase of the states
syms A_1 A_2 A_3 A_4 A_9 A_10 A_11 real;
syms phi_1 phi_2 phi_3 phi_4 phi_9 phi_10 phi_11 real;
Amplitudes = [A_1; A_2; A_3; A_4; A_9; A_10; A_11; A_alpha; w*A_alpha; w*H*b];
phis = [phi_1; phi_2; phi_3; phi_4; phi_9; phi_10; phi_11; 0; pi/2; phi+pi/2];

% Averaged increases/decreases
dCn = averaged_Taylor(Cn,Q,Amplitudes,phis);
% dCd = averaged_Taylor(Cd,Q,Amplitudes,phis); toc 