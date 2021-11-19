% University of California, Irvine - Summer 2021
% Laura Pla Olea - lplaolea@uci.edu

clc; close all; clear;

%% Inputs

n = 2; % Pullback order
nT = 2; % Order of the Taylor expansion
norder = 2; % Order of the averaged output that we want to keep

In = 'HA'; % Motion: 'H' for plunging, 'A' for pitching, 'HA' for both
Out = 'L'; % Output: 'L' for lift coefficient, 'D' for drag coefficient

%% Definitions

% States of the system
x = sym('x',[11 1]); % Original states of the BL system
x(5:8) = [];
syms a da dh real; % Angle of attack, pitching velocity, plunging velocity

% Parameters of the system
syms t real; % Time - Independent variable
syms rho real; % Density
syms beta M U b real; % Compressibility factor, Mach number, velocity, airfoil semichord
syms A1 A2 b1 b2; % Circulatory component constants
syms Ka Kq Ti real; % Noncirculatory component constants
syms Tp Tf Tv real; % Nonlinear flow components
syms fp; % Effective point of separation
syms w A_alpha H phi real; % Angular velocity, pitching amplitude, plunging amplitude, phase between pitching and plunging

syms CNa;
syms CNs(a_eff_); % Steady lift coefficient (function of alpha and plunging)
assume(CNs(a_eff_),'real');

syms dCv(x9,x10,aE_,a_eff_,da_); % Derivative of the vortex strength
assume(dCv(x9,x10,aE_,a_eff_,da_),'real');

syms eps real; % Epsilon (to determine the order of the terms)

% Order of the variables
w = w/eps;
U = U/eps;
A_alpha = A_alpha*eps;
H = H*eps;
if strcmp(In,'H')
    A_alpha = 0; phi = 0;
elseif strcmp(In,'A')
    H = 0;
end

c = 2*b; % Airfoil chord
a_eff = a+atan(dh/U); % Effective angle of attack
aE = (2*U/c)*beta^2*(A1*b1*x(1)+A2*b2*x(2)); % Effective angle of attack
fp = (2*sqrt(CNs(x(5)/CNa)/x(5))-1)^2; % Effective point of separation

%% State Space Model

% States of the system
Q = [x; a; da; dh];

% Drift vector
f = [-b1*beta^2*(2*U/c)*x(1)+a_eff+da*c/(2*U);
    -b2*beta^2*(2*U/c)*x(2)+a_eff+da*c/(2*U);
    -x(3)/(Ka*Ti)+a_eff;
    -x(4)/(Kq*Ti)+da*c/U;
    (CNa*beta^2*(2*U/c)*(A1*b1*x(1)+A2*b2*x(2))-4/(M*Ka*Ti)*x(3)-1/(M*Kq*Ti)*x(4)-x(5)+4/M*a_eff+da*c/U/M)*2*U/c/Tp;
    (-x(6)/Tf+fp/Tf)*2*U/c;
    (-x(7)/Tv+dCv(x(5),x(6),aE,a_eff,da)/Tv)*2*U/c;
    da;
    0;
    0];
f = formula(f);

ga = [0 0 0 0 0 0 0 0 1 0].'; % Control vector of the pitching motion
gh = [0 0 0 0 0 0 0 0 0 1].'; % Control vector of the plunging motion

%% Lift and drag coefficients

q = da*c/U; % Pitch rate

% Normal force coefficients
Cnia = -4/M*x(3)/(Ka*Ti)+4/M*a_eff; % Non-circulatory load due to the angle of attack
Cniq = -1/M*x(4)/(Kq*Ti)+q/M; % Non-circulatory load due to the pitch rate
Cni = Cnia+Cniq; % Total non-circulatory (impulsive) normal force coefficient

Cnf = CNa*aE*((1+sqrt(x(6)))/2)^2; % Circulatory load due to TE separation
Cnv = x(7); % Circulatory load due to LE separation (vortex)

Cn = Cnf+Cnv+Cni; % Normal force coefficient

% Chord force coefficient
syms eta real; % Parameter that accounts for the viscous effects
Ccp = CNa*aE^2; % Chord force under attached flow conditions
Cc = eta*Ccp*sqrt(x(6)); % Chord force coefficient under TE separation

% Total lift and drag coefficients
Cl = Cn*cos(a)+Cc*sin(a); % Lift coefficient
Cd = Cn*sin(a)-Cc*cos(a); % Drag coefficient

%% Pullback

% Inputs
ua = w^2*A_alpha*cos(w*t);
uh = w^2*H*b*cos(w*t+phi);

G = ga*ua+gh*uh;

% Pullback of f along the flow of G
F = simplify(pullback(Q,f,G,t,n),'IgnoreAnalyticConstraints',true);
Gvec_avg = simplify(w/(2*pi)*int(F-f,t,0,2*pi/w));
F_avg = simplify(f+Gvec_avg);

%% Equilibrium of the averaged dynamics

dheq = 0; % The equilibrium of dh is automatically satisfied, so for simplification we assume it's zero

% Solving the equilibrium points of the averaged dynamics
daeq = solve(F_avg(8)==0,da); % da of equilibrium
Favg_sol = subs(F_avg,[da dh],[daeq dheq]);

xeq = x;
for i = 1:length(x)
    xeq(i) = solve(Favg_sol(i)==0,x(i));
    if i==5
        xeq(i) = simplify(subs(xeq(i),A2,1-A1));
    end
    Favg_sol = subs(Favg_sol,x(i),xeq(i));
end

Qeq = [xeq; a; daeq; dheq];

%% Averaged output

% Amplitude and phase of the states
Ax = sym('Ax',[11 1]); % Amplitude of the original BL states
Ax(5:8) = [];
phix = sym('phix',[11 1]); % Phase of the original BL states
phix(5:8) = [];

dAx = Ax.*cos(w*t+phix);
ddalpha = w*A_alpha*sin(w*t);
dalpha = -A_alpha*cos(w*t);
dhdot = w*H*b*(sin(w*t+phi)-sin(phi));
dq = [dAx; dalpha; ddalpha; dhdot];

if strcmp(Out,'L')
    Psi = Cl;
elseif strcmp(Out,'D')
    Psi = Cd;
end

% Taylor expansion
Psi_Taylor = Taylor_expansion(Psi,Q,Qeq,dq,2);
Psi_avg = simplify(int(Psi_Taylor,t,0,2*pi/w)*w/(2*pi));

%% Order of terms

epsmax = polynomialDegree(expand(Psi_avg),eps);
ord_terms = simplify(coeffs(formula(expand(Psi_avg)),eps,'All'));

for i = 1:epsmax+1
    disp(['O(',char(eps^(i-1)),')']);
    pretty(ord_terms(epsmax+2-i));
end

if norder>epsmax
    norder = epsmax;
end
Psi_truncated = 0;
for i = 0:norder
    Psi_truncated = Psi_truncated+ord_terms(epsmax+1-i);
end