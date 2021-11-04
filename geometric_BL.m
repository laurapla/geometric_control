% University of California, Irvine - Summer 2021
% Laura Pla Olea - lplaolea@uci.edu
tic
clc; close all; clear;

%% State Space

syms x1 x2 x3 x4 x9 x10 x11 real; % Original states of the system
syms a da dh real; % Angle of attack, pitching velocity, plunging velocity
Q = [x1 x2 x3 x4 x9 x10 x11 a da dh].'; % States of the system

syms CNa;
% syms CNs(x1,x2,x3,x4,x9,x10,a,da,dh); % Static normal coefficient
% assume(CNs(x1,x2,x3,x4,x9,x10,a,da,dh),'real');
% syms CNa(x1,x2,x3,x4,x9,x10,a,da,dh); % Slope of the static normal coefficient
% assume(CNa(x1,x2,x3,x4,x9,x10,a,da,dh),'real');
syms beta M V b real; % Compressibility factor, Mach number, velocity, airfoil semichord
syms A1 A2 b1 b2; % Circulatory component constants
syms Ka Kq Ti real; % Noncirculatory component constants
syms Tp Tf Tv; % Nonlinear flow components
syms fp(x1,x2,x3,x4,x9,x10,a,da,dh);
assume(fp(x1,x2,x3,x4,x9,x10,a,da,dh),'real');
syms dCv(x1,x2,x3,x4,x9,x10,a,da,dh);
assume(dCv(x1,x2,x3,x4,x9,x10,a,da,dh),'real');

% fp = (2*sqrt(CNs/CNa)-1)^2; % Effective point of separation
c = 2*b; % Airfoil chord

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

%% Pullback

n = 2; % Pullback order

syms w A_alpha H phi real; % Constants of the inputs
syms t real; % Independent variable

% Inputs
ua = w^2*A_alpha*cos(w*t);
uh = w^2*H*b*cos(w*t+phi);
% G = ga*ua+gh*uh;
G = [ga gh]; Us = [ua uh];

% Pullback of f along the flow of G
F = simplify(pullback(Q,f,G,Us,t,n),'IgnoreAnalyticConstraints',true);
Gvector_averaged = simplify(w/(2*pi)*int(F-f,t,0,2*pi/w));


% %% Series expansion (averaged values)
% 
% % Amplitude and phase of the states
% syms A_1 A_2 A_3 A_4 A_9 A_10 A_11 real;
% syms phi_1 phi_2 phi_3 phi_4 phi_9 phi_10 phi_11 real;
% Amplitudes = [A_1; A_2; A_3; A_4; A_9; A_10; A_11; A_alpha; w*A_alpha; w*H*b];
% phis = [phi_1; phi_2; phi_3; phi_4; phi_9; phi_10; phi_11; 0; pi/2; phi+pi/2];
% 
% % Averaged increases/decreases
% dCn = averaged_Taylor(Cn,Q,Amplitudes,phis);
% % dCd = averaged_Taylor(Cd,Q,Amplitudes,phis);
toc 