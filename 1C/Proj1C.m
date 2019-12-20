%% Progetto 1C - Link flessibile
%%
% Corradini, di Nuzzo, Frick, Ragazzini, Zappacosta
% Gruppo A
%

format compact;

%% Descrizione e requisiti del sistema
% Nella nuova applicazione della azienda che commissiona il progetto 
% si prevede di utilizzare una struttura meccanica particolarmente leggera. 
% Questa però presenta il problema di una flessibilità non trascurabile 
% intrinseca nei componenti meccanici utilizzati. 
% Ciò rende difficile il posizionamento dell’estremità non attuata in una 
% posizione fissa desiderata.
%
% Per l’applicazione che l’azienda ha in mente si devono rispettare per il 
% sistema linearizzato determinate caratteristiche:
%
% # Errore a regime nullo con riferimento a gradino con ampiezza $w(t) = W sca(t)$.
% # Per garantire una certa robustezza del sistema si deve avere un margine di fase $\phi_m \geq 45^\circ$.
% # Il sistema può accettare una sovraelongazione percentuale al massimo dell’1\% : $S\% \leq 1\%$.
% # Tempo di assestamento all'1\% $T_{a1} = 0.8$ (opzionalmente 0.4).
% # Abbattimento del rumore di 20 volte.
%
%
% Il rumore si fa sentire a $\omega_n > 200 rad/s$ con ampiezza $A_n = 0.05$.

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
% Il requisito sulla sovraelongazione si traduce in un requisito sul
% margine di fase che è maggiore di quello specificato, quindi si prenderà 
% in considerazione solo questo.
xi = sqrt(log(S_max)^2/(pi^2+log(S_max)^2));
phi_m = xi * 100

%% Definizione del regolatore
%
% # Per avere errore a regime nullo è necessario che $L(s)$ abbia un polo
% nell'origine, ma $G(s)$ ne ha due, che inoltre abbassano di molto la fase: 
% si progetta quindi $R(s)$ in modo che abbia uno zero vicino all'origine
% che cancelli il polo;
% # Ci si serve di due reti anticipatrici, una
% con punto medio in $\omega_1 = 1/ (T_1 \sqrt{\alpha_1}) \approx 5 \cdot 10^{-1}$
% e un'altra con punto medio in $\omega_2 = 1/ (T_2 \sqrt{\alpha_2})
% \approx 7 \cdot 10^{-3}$;
% # La prima rete anticipatrice ha una larghezza di banda di circa una 
% decade, uno zero a $\omega_{z1} \approx -0.06$ e un polo a $omega_{p1}
% \approx -70$.
% # La seconda rete anticipatrice ha una larghezza di banda di circa tre
% decadi, uno zero a $\omega_{z2} \approx -20$ e un polo a $omega_{p2}
% \approx -10^{-3}$.
%

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
% Il sistema rispetta sia le specifiche sulla sovraelongazione sia quelle
% sul tempo di assestamento all'1%.
stepinfo(F, 'SettlingTimeThreshold',0.01)

%% Simulink
% Il sistema non linearizzato è stabile per valori dell'ingresso $w$ minori
% di $W/8 sca(t)$
w_lin = W;
lin_sim = sim('Simul1C');
w_nonlin = W/8;
nonlin_sim = sim('NonLin1C');

y_lin = lin_sim.get('y');
y_nonlin = nonlin_sim.get('y');

figure;
hold on;
xlim([0 100]);
plot (y_lin);
title("Linearized system - step response with w = W");
hold off;

figure;
hold on;
xlim([0 200]);
plot (y_nonlin);
title("Non linearized system - step response with w = W/8");
hold off;

