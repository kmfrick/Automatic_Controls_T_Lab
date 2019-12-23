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
B_n_dB = 20 * log10(B_n);
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
G_mat=tf(N,D);
G = 1000 * (s+15) / (s^2 * (s^2 + 30 * s + 450));

%% Definizione dell'intervallo di frequenze del diagramma di Bode
omega_plot_min=10^(-2);
omega_plot_max=10^5;

%% Requisiti sul margine di fase
xi = sqrt(log(S_max)^2/(pi^2+log(S_max)^2));
phi_m = xi * 100;
omega_c_opt = 460 / (phi_m * T_a0);
omega_c_stand = 460 / (phi_m * T_a1);

%% Definizione del regolatore per specifiche standard
% Zero nell'origine per alzare la fase a valori accettabili
margin(s * G); grid on;
% Si ha margine di fase > 83 gradi per omega < 10
gain_sta = 10^(-log10(abs(evalfr(s * G, 9.5 * 1i))));
% Si aggiunge un polo lontano dall'origine per la fisica realizzabilità
R_sta = gain_sta * s / (1 + s * 1e-5);
L_sta = R_sta * G;
F_sta = L_sta / (1 + L_sta);


%% Definizione del regolatore ottimale
gain = 3.4e-2;
T_lead_1 = 17.42;
alpha_lead_1 = 5.741e-5;
T_lead_2 = 0.05;
alpha_lead_2 = 0.286;
R_lead_1 = (1  + T_lead_1 * s) / (1 + alpha_lead_1 * T_lead_1 * s);
R_lead_2 = (1  + T_lead_2 * s) / (1 + alpha_lead_2 * T_lead_2 * s);
R_opt = gain * R_lead_1 * R_lead_2;
L_opt = R_opt * G;
F_opt = L_opt / (1 + L_opt);

%% Diagramma di Bode del sistema in anello aperto
figure;
hold on;
patch([omega_plot_min,omega_c_opt,omega_c_opt, omega_plot_min],[0, 0, -300, -300],'green','FaceAlpha',0.3,'EdgeAlpha',0);
patch([omega_plot_min,omega_c_stand,omega_c_stand, omega_plot_min],[0, 0, -300, -300],'green','FaceAlpha',0.5,'EdgeAlpha',0);
patch([omega_n,omega_plot_max,omega_plot_max,omega_n],[-B_n_dB,-B_n_dB,100,100],'red','FaceAlpha',0.5,'EdgeAlpha',0);
[Mag, phase, omega] = bode(G, {omega_plot_min, omega_plot_max});
margin(Mag, phase, omega); grid on;
patch([omega_c_stand,omega_n, omega_n, omega_c_stand], [-180+phi_m,-180+phi_m, -360, -360], 'blue','FaceAlpha',0.5,'EdgeAlpha',0);

%% Diagramma di Bode del sistema con regolatore standard
figure;
hold on;
patch([omega_plot_min,omega_c_opt,omega_c_opt, omega_plot_min],[0, 0, -200, -200],'green','FaceAlpha',0.3,'EdgeAlpha',0);
patch([omega_plot_min,omega_c_stand,omega_c_stand, omega_plot_min],[0, 0, -200, -200],'green','FaceAlpha',0.5,'EdgeAlpha',0);
patch([omega_n,omega_plot_max,omega_plot_max,omega_n],[-B_n_dB,-B_n_dB,100,100],'red','FaceAlpha',0.5,'EdgeAlpha',0);
[Mag, phase, omega] = bode(L_sta, {omega_plot_min, omega_plot_max});
margin(Mag, phase, omega);  
grid on;
patch([omega_c_stand,omega_n, omega_n, omega_c_stand], [-180+phi_m,-180+phi_m, -360, -360], 'blue','FaceAlpha',0.5,'EdgeAlpha',0);


%% Diagramma di Bode del sistema con regolatore ottimale
figure;
hold on;
patch([omega_plot_min,omega_c_opt,omega_c_opt, omega_plot_min],[0, 0, -200, -200],'green','FaceAlpha',0.3,'EdgeAlpha',0);
patch([omega_plot_min,omega_c_stand,omega_c_stand, omega_plot_min],[0, 0, -200, -200],'green','FaceAlpha',0.5,'EdgeAlpha',0);
patch([omega_n,omega_plot_max,omega_plot_max,omega_n],[-B_n_dB,-B_n_dB,100,100],'red','FaceAlpha',0.5,'EdgeAlpha',0);
[Mag, phase, omega] = bode(L_opt, {omega_plot_min, omega_plot_max});
margin(Mag, phase, omega);  
grid on;
patch([omega_c_stand,omega_n, omega_n, omega_c_stand], [-180+phi_m,-180+phi_m, -360, -360], 'blue','FaceAlpha',0.5,'EdgeAlpha',0);

%% Risposta allo scalino del sistema in anello chiuso con regolatore standard
figure;
hold on;
step(F_sta);
stepinfo(F_sta, 'SettlingTimeThreshold',0.01)

%% Risposta allo scalino del sistema in anello chiuso con regolatore ottimale
figure;
hold on;
step(F_opt);
stepinfo(F_opt, 'SettlingTimeThreshold',0.01)

%% Simulink
R_sim = R_opt;
w_lin = W;
w_nonlin = W/8;
lin_sim = sim('Simul1C');
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

