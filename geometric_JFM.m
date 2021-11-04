% University of California, Irvine - Winter 2020
% Laura Pla Olea - lplaolea@uci.edu

clc; close all; clear;

%% Inputs

nc = 2; % Number of circulatory states
n = 2; % Pullback order

%% Definitions

% States of the system
syms xv a da dh real; % Viscous state, angle of attack, derivative of a, plunging speed
xc = sym('xc', [nc 1]); % Circulatory states

% Parameters of the system
syms t real; % Time - Independent variable
syms U b e real; % Forward speed, airfoil semichord, pitching axis
syms mv tv kda real; % Added mass, added mass constant, da constant
syms w A_alpha H phi real; % Angular velocity, pitching amplitude, plunging amplitude, phase between pitching and plunging

syms CLs(a_eff); % Steady lift coefficient (function of alpha and plunging)
assume(CLs(a_eff),'real');

syms eps; % Epsilon (to determine the order of the terms)

% Transfer function between the quasi-steady circulation and the effective
% circulation
syms kdc khf; % DC gain, high-frequency gain
aa = sym('a', [n 1]);
bb = sym('b', [n 1]);
bb(1) = kdc*aa(1);
A = [zeros(1,nc-1) -aa(1);
    eye(n-1) -aa(2:n)];
B = bb-khf*aa;
C = [zeros(1,n-1) 1];

c = 2*b;
ka = U*c/2;
a_eff = a+atan(dh/U);
Gamma0 = ka*CLs(a_eff)+kda*da; % Quasi-steady circulation

%% State Space Model

% States of the system
Q = [xc; xv; a; da; dh].';

% Drift
f = [A*xc+B*Gamma0; -xv/tv+(U*cos(a)-dh*sin(a))*da/tv; da; 0; 0];
f = formula(f);

% Control vectors
ga = [zeros(1,n) -e/tv 0 1 0].'; % Angle of attack
gh = [zeros(1,n) cos(a)/tv 0 0 1].'; % Plunging

%% Pullback

% Inputs
ua = w^2*A_alpha*cos(w*t);
uh = w^2*H*b*cos(w*t+phi);

G = ga*ua+gh*uh;

% Pullback of f along the flow of G
F = simplify(pullback(Q,f,G,t,n),'IgnoreAnalyticConstraints',true);
Gvector_averaged = simplify(w/(2*pi)*int(F-f,t,0,2*pi/w));