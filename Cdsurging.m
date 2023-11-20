% University of California, Irvine - Fall 2023
% Laura Pla Olea - lplaolea@uci.edu

clc; clear; close all;

%% Get data

% Static
static_data = 'C:\Users\laura\Documents\UCI\JFM Paper\CL figure\Results_3.xlsx';
alpha_sim = readmatrix(static_data,'Sheet','Sheet1','Range','A2:A30');
Cls_sim = readmatrix(static_data,'Sheet','Sheet1','Range','B2:B30');

ks_data = 'C:\Users\laura\Documents\GitHub\geometric_control\ks-alpha.txt';
ks_data_mat = readmatrix(ks_data);
ks_alpha = ks_data_mat(:,1);
ks_dat = ks_data_mat(:,2);

%% Parameters

khf = 1/2;
b = 1;
tau_v = 2e2;
Re = 5e5;
rho = 1.225;
mu = 1.789e-5;

%% Static lift coefficient

xx = linspace(min(alpha_sim),max(ks_alpha),1e3);
yy = spline(alpha_sim,Cls_sim,xx);
yyy = spline(ks_alpha,ks_dat,xx);
% plot(alpha,Cls,'o',xx,yy)

alpha = xx;
Cls = yy;
ks = 2*pi*yyy;

%% Cd terms

U = Re*mu/(rho*2*b);

ord1 = (1-khf)*(Cls.*tand(alpha)-ks/(2*pi^2).*Cls.^2-pi*b/(U*tau_v)*sin(alpha).^2);
ord2 = (1-khf)*(Cls.*tand(alpha)-ks/(2*pi^2)*3*(7-3*khf)/4.*Cls.^2);

%% Figure

line_width = 1.7;
font_lgd = 14;
font_labels = 14;
inter = 2;

symbols = {'-', '--', ':','-.','-', '--', ':','-.','-', '--', ':','-.','-', '--', ':','-.','-', '--', ':','-.'};
Okabe_Ito = [0.902 0.624 0; 0.337 0.737 0.914; 0 0.62 0.451;
    0.941 0.894 0.259; 0 0.447 0.698; 0.835 0.369 0; 0.8 0.475 0.655];

figure;
colororder(Okabe_Ito)

yyaxis left
plot(alpha,Cls,'-','LineWidth',line_width);
ylabel('$C_{L,s}$','interpreter','latex','FontSize',font_labels);

yyaxis right
plot(alpha,ord1,'--',alpha,ord2,':','LineWidth',line_width)
xlabel('$\alpha^{*}, ^{\circ}$','interpreter','latex','FontSize',font_labels);
ylabel('$\chi$','interpreter','latex','FontSize',font_labels);
legend('','$\chi_{\sigma}$','$\chi_{\sigma^{2}}$','Location','best','interpreter','latex','FontSize',font_lgd)
grid on;