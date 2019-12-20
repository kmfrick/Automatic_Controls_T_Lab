%% Progetto 1C - Link flessibile
%%
% Corradini, di Nuzzo, Frick, Ragazzini, Zappacosta
% Gruppo A
%

format compact;

%% Definizione dei parametri del sistema
g = 9.81;
K = 3;
ro = 0.2;
L = 1;
I = 0.01;
J = 0.02;
M = 0.06;
h = 1;
T_a1 = 0.8;
T_a0 = 0.4;
W = 0.0872;
x_1_ref = pi/2;
x_2_ref = 0;
x_3_ref = M * g * L / K + x_1_ref;
x_4_ref = 0;
u_ref = M * g * L;
S_max = 0.01;
omega_n = 200;
A_n = 0.05 * pi / 180;
B_n = 20;
B_n_dB = 20 * log(20) / log(10);
y_ref = pi/2;

%% Matrici linearizzate
A = [0      1       0       0;
    (-K/J)  (-ro/J) K/J   ro/J;
    0       0       0       1;
    K/I   ro/I  (-K/I) (-ro/I)];
B = [0; 0; 0; 1/I];
C = [1 0 0 0];
D = 0;

%% Definizione funzione di trasferimento
s=tf('s');
[N,D]=ss2tf(A,B,C,D);
G=tf(N,D);

%% Definizione dell'intervallo di frequenze del diagramma di Bode
omega_plot_min=10^(-2);
omega_plot_max=10^5;

%% Requisiti sul margine di fase
xi = sqrt(log(S_max)^2/(pi^2+log(S_max)^2));
phi_m = xi * 100;

%% Definizione del regolatore
gain = 3.4e-2;
R_s = gain;
T_lead_1 = 17.42;
alpha_lead_1 = 8.2e-4;
T_lead_2 = 0.05;
alpha_lead_2 = 0.02;
R_lead_1 = (1  + T_lead_1 * s) / (1 + alpha_lead_1 * T_lead_1 * s);
R_lead_2 = (1  + T_lead_2 * s) / (1 + alpha_lead_2 * T_lead_2 * s);
R = R_s * R_lead_1 * R_lead_2;

L_ = R * G;
F = L_ / (1 + L_);

%% Diagramma di Bode del sistema in anello aperto
figure;
patch([omega_n,omega_plot_max,omega_plot_max,omega_n],[-B_n_dB,-B_n_dB,100,100],'red','FaceAlpha',0.3,'EdgeAlpha',0);
hold on; 
[Mag, phase, omega] = bode(G, {omega_plot_min, omega_plot_max});
margin(Mag, phase, omega); grid on;

%% Diagramma di Bode del sistema con regolatore
figure;
hold on;
patch([omega_n,omega_plot_max,omega_plot_max,omega_n],[-B_n_dB,-B_n_dB,100,100],'red','FaceAlpha',0.3,'EdgeAlpha',0);
[Mag, phase, omega] = bode(L_, {omega_plot_min, omega_plot_max});
margin(Mag, phase, omega);  
grid on;

%% Risposta allo scalino del sistema in anello chiuso con regolatore
stepinfo(F, 'SettlingTimeThreshold',0.01)

%% Simulink
w_lin = W;
lin_sim = sim('Simul1C');
w_nonlin = W/8;
nonlin_sim = sim('NonLin1C');

y_lin = lin_sim.get('y');
y_nonlin = nonlin_sim.get('y');

figure;
hold on;
xlim([0 4]);
plot (y_lin);
title("Linearized system - step response with w = W");
hold off;

figure;
hold on;
xlim([0 4]);
plot (y_nonlin);
title("Non linearized system - step response with w = W/8");
hold off;

figure;
hold on;
xlim([0 80]);
plot (y_nonlin);
title("Non linearized system - step response with w = W/8");
hold off;

