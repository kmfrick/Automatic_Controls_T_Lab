%% Project 1C
% Corradini, di Nuzzo, Frick, Ragazzini, Zappacosta
% Group A

% x_dot=f(x,u,t)

% x_dot_1 = x_2;
% x_dot_2 = (-M g L sin (x1) - K (x_1 - x_3) - ro (x_2 - x_4)) / J
% x_dot_3 = x_4
% x_dot_4 = (K(x_1 - x_2) + ro (x_2 - x_4) + u) / I

% y = x_1

%% Project specs:
% zero steady state error with a step reference signal of 5
% Mf > 45 deg
% S_% < 1%
% T_a1 = 0.8 (0.4 optionally)
% Noise reduction of about 20 times
% Noise characteristics: omega_n > 200 [rad/s], amplitude=0.05

omega_n = 200;
A_n = 0.05;
B_n = 20;
y_ref = pi/2;

% System parameters definition

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
W = 5;
x_1_ref = pi/2;
x_2_ref = 0;
x_3_ref = M * g * L / K + x_1_ref;
x_4_ref = 0;
u_ref = M * g * L;
S_max = 0.01;

% Linearized matrices
A = [0      1       0       0;
    (-K/J)  (-ro/J) K/J   ro/J;
    0       0       0       1;
    K/I   ro/I  (-K/I) (-ro/I)];
B = [0; 0; 0; 1/I];
C = [1 0 0 0];
D = 0;

% Transfer function definition

s=tf('s');
[N,D]=ss2tf(A,B,C,D);
G=tf(N,D);

% Definition of Bode diagram frequency range
omega_plot_min=10^(-2);
omega_plot_max=10^5;

% Plot Bode diagram of open-loop system
figure;
patch([omega_n,omega_plot_max,omega_plot_max,omega_n],[-B_n,-B_n,100,100],'red','FaceAlpha',0.3,'EdgeAlpha',0);
hold on; 
[Mag, phase, omega] = bode(G, {omega_plot_min, omega_plot_max});
margin(Mag, phase, omega); grid on;

xi = sqrt(log(S_max)^2/(pi^2+log(S_max)^2));
% The overshoot request maps into a phase margin request:
phi_m=xi*100;
% Mf>83% this is a stronger request that the robust one (Mf>45%),
% so we will only take this into account.

% Define regulator
% Zero steady state error: requires that L has a pole in the origin
alpha_lead = 0.5;
T_lead = 3e-2;
gain = 3.3e-1;
R_s = s * gain /  (1 + 1e-6 * s);
R_lead = (1  + T_lead * s) / (1 + alpha_lead * T_lead * s);
R_d = R_lead;
R = R_s * R_d;
L = R * G;
figure;
hold on;
patch([omega_n,omega_plot_max,omega_plot_max,omega_n],[-B_n,-B_n,100,100],'red','FaceAlpha',0.3,'EdgeAlpha',0);
[Mag, phase, omega] = bode(L, {omega_plot_min, omega_plot_max});
margin(Mag, phase, omega);  
grid on;

% Plot step response of regulated closed-loop system

F = L / (1 + L);
figure;
hold on;
step(F);
grid on;

open('Simul1C.slx')