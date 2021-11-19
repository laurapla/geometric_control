% University of California, Irvine - Fall 2021
% Laura Pla Olea - lplaolea@uci.edu

clc; close all; clear;

%% Inputs

nc = 4; % Number of circulatory states
n = 2; % Pullback order
nT = 2; % Order of the Taylor expansion
norder = 2; % Order of the averaged output that we want to keep

In = 'HA'; % Motion: 'H' for plunging, 'A' for pitching, 'HA' for both
Out = 'L'; % Output: 'L' for lift coefficient, 'D' for drag coefficient, 'S' for point of separation

%% Definitions

% States of the system
syms xv a da dh real; % Viscous state, angle of attack, derivative of a, plunging speed
xc = sym('xc', [nc 1]); % Circulatory states

% Parameters of the system
syms t real; % Time - Independent variable
syms rho real; % Density
syms U b e real; % Forward speed, airfoil semichord, pitching axis
syms mv tv kda real; % Added mass, added mass constant, da constant
syms w A_alpha H phi real; % Angular velocity, pitching amplitude, plunging amplitude, phase between pitching and plunging

syms CLs(a_eff); % Steady lift coefficient (function of alpha and plunging)
assume(CLs(a_eff),'real');

syms eps real; % Epsilon (to determine the order of the terms)

% Transfer function between the quasi-steady circulation and the effective
% circulation
syms kdc khf real; % DC gain, high-frequency gain
aa = sym('a', [nc 1]);
bb = sym('b', [nc 1]);
bb(1) = kdc*aa(1);
A = [zeros(1,nc-1) -aa(1);
    eye(nc-1) -aa(2:nc)];
B = bb-khf*aa;
C = [zeros(1,nc-1) 1];

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

c = 2*b;
ka = U*c/2;
a_eff = a+atan(dh/U);
Gamma0 = ka*CLs(a_eff)+kda*da; % Quasi-steady circulation

%% State Space Model

% States of the system
Q = [xc; xv; a; da; dh];

% Drift
f = [A*xc+B*Gamma0; -xv/tv+(U*cos(a)-dh*sin(a))*da/tv; da; 0; 0];
f = formula(f);

% Control vectors
ga = [zeros(1,nc) -e/tv 0 1 0].'; % Angle of attack
gh = [zeros(1,nc) cos(a)/tv 0 0 1].'; % Plunging

if strcmp(Out,'S')
    
    % States of the system
    syms xs;
    Q = [Q(1:nc+1); xs; Q(nc+2:end)];
    
    % Parameters
    syms tau1 tau2; % Point of separation, point of separation constants
    x0 = (sqrt(2*CLs(a_eff-tau2*da)/(pi*(a_eff-tau2*da)))-1)^2;
    
    % Drift
    f = [f(1:nc+1); -xs/tau1+x0/tau1; f(nc+2:end)];
    
    % Control vectors
    ga = [ga(1:nc+1); 0; ga(nc+2:end)];
    gh = [gh(1:nc+1); 0; gh(nc+2:end)];
    
end

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
if strcmp(Out,'S')
    daeq = solve(F_avg(nc+3)==0,da); % da of equilibrium
    Favg_sol = subs(F_avg,[da dh],[daeq dheq]);
    xseq = solve(Favg_sol(nc+2)==0,xs); % xs of equilibrium
else
    daeq = solve(F_avg(nc+2)==0,da); % da of equilibrium
    Favg_sol = subs(F_avg,[da dh],[daeq dheq]);
end
xveq = solve(Favg_sol(nc+1)==0,xv); % xv of equilibrium
xceq = -inv(A)*(Favg_sol(1:nc)-A*xc); % xc of equilibrium

Qeq = [xceq; xveq; a; daeq; dheq];

if strcmp(Out,'S')
    Qeq = [Qeq(1:nc+1); xseq; Qeq(nc+2:end)];
end

%% Averaged output

Ax = sym('Ax',[nc 1]); % Amplitude of the circulatory states
Ax = Ax*eps;
phix = sym('phix',[nc 1]); % Phase of the circulatory states
syms Av phiv real; % Amplitude and phase of the viscous state
Av = Av*eps;

dAx = Ax.*cos(w*t+phix);
dAv = Av*cos(w*t+phiv);
ddalpha = w*A_alpha*sin(w*t);
dalpha = -A_alpha*cos(w*t);
dhdot = w*H*b*(sin(w*t+phi)-sin(phi));
dq = [dAx; dAv; dalpha; ddalpha; dhdot];

if ~strcmp(Out,'S')
    
    % Lift force
    L = rho*U*(C*xc+khf*Gamma0)+mv*xv*cos(a);
    
    if strcmp(Out,'L')
        Psi = L/(rho*U^2*b);
    elseif strcmp(Out,'D')
        syms ks(a_eff)
        a_eff = a+atan(dh/U);
        ks = ks(a_eff);
        % Drag force
        D = L*tan(a)-ks*rho*b*((C*xc+khf*Gamma0)/(2*pi*b)-b*da/2)^2;
        Psi = D/(rho*U^2*b);
    end
    
    % Taylor expansion
    Psi_Taylor = Taylor_expansion(Psi,Q,Qeq,dq,2);
    Psi_avg = simplify(int(Psi_Taylor,t,0,2*pi/w)*w/(2*pi));

else
    
    Psi_avg = xseq;
    
end

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