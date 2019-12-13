%% Progetto 1C - Link flessibile
%%
% Corradini, di Nuzzo, Frick, Ragazzini, Zappacosta
% Gruppo A

% x_dot=f(x,u,t)

% x_dot_1 = x_2;
% x_dot_2 = (-M g L sin (x1) - K (x_1 - x_3) - ro (x_2 - x_4)) / J
% x_dot_3 = x_4
% x_dot_4 = (K(x_1 - x_2) + ro (x_2 - x_4) + u) / I

% y = x_1

format compact;

%% Specifiche
%
% # Errore a regime nullo con ingresso $w = A sca(t)$
% # Margine di fase $\phi_m > 45^{\circ}$
% # Sovraelongazione percentuale massima $S_\% < 1 \%$
% # Tempo di assestamento all'1\% $T_{a1} = 0.8$ (opzionalmente $0.4$)
% # Abbattimento del rumore di 20 volte
% # Caratteristiche del rumore: $\omega_n > 200 rad/s, A_n = 0.05$
%

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
A_n = 0.05;
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
zpk(G)

%% Definizione dell'intervallo di frequenze del diagramma di Bode
omega_plot_min=10^(-2);
omega_plot_max=10^5;

%% Diagramma di Bode del sistema in anello aperto
figure;
patch([omega_n,omega_plot_max,omega_plot_max,omega_n],[-B_n_dB,-B_n_dB,100,100],'red','FaceAlpha',0.3,'EdgeAlpha',0);
hold on; 
[Mag, phase, omega] = bode(G, {omega_plot_min, omega_plot_max});
margin(Mag, phase, omega); grid on;

%% Requisiti sul margine di fase
% Il requisito sulla sovraelongazione si traduce in un requisito sul
% margine di fase che è maggiore di quello specificato, quindi si prenderà 
% in considerazione solo questo.
xi = sqrt(log(S_max)^2/(pi^2+log(S_max)^2));
phi_m = xi * 100


%% Definizione del regolatore
%
% # Per avere errore a regime nullo è necessario che $L(s)$ abbia un polo
% nell'origine, ma $G(s)$ ne ha due che abbassano di molto la fase: 
% si progetta quindi $R(s)$ in modo che abbia uno zero nell'origine.
% # La fase è molto negativa, quindi si usa una rete anticipatrice
% $R_{lead}$ per alzare la fase a $\omega \approx 12.9$.
% # Si inserisce poi un polo a sei decadi dall'origine per garantire la
% fisica realizzabilità.
%
alpha_lead = 0.2;   
T_lead = 3e-2;
gain = 2.4e-3;
R_s = gain;
R_lead = (1  + T_lead * s) / (1 + alpha_lead * T_lead * s);
R_d = (s * 1e2 + 1) * R_lead /  (1 + 1e-6 * s);
R = R_s * R_d;
L_ = R * G;
zpk(R)
zpk(L_)
figure;
hold on;
patch([omega_n,omega_plot_max,omega_plot_max,omega_n],[-B_n_dB,-B_n_dB,100,100],'red','FaceAlpha',0.3,'EdgeAlpha',0);
[Mag, phase, omega] = bode(L_, {omega_plot_min, omega_plot_max});
margin(Mag, phase, omega);  
grid on;

%% Risposta in frequenza del sistema in anello chiuso con regolatore
% Il sistema rispetta sia le specifiche sulla sovraelongazione sia quelle
% sul tempo di assestamento all'1%.
F = L_ / (1 + L_);
figure;
hold on;
step(F);
figure;
[Mag, phase, omega] = bode(L_, {omega_plot_min, omega_plot_max});
margin(Mag, phase, omega);  
grid on;

%open('Simul1C.slx')