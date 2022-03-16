% Cankan ÇELÝK 2226736 ME442 Term PROJECT
%% System Properties
clc; clear; close all;
%System Properties
N=100; %Turns
L=10; %H
R=500; %Ohm
Io=10; %A
m=100; %kg
A=20; %m^2
nu=2.9*10^-4; %F/m
g = 9.81; % m^2/s
K=(nu*A*(N^2))/2;  %Linearized Constant
Zo=1.72; %m
inv=-1; %AC->DC
%###########
%Block Representation
b1_num_A=[inv]; b1_den_A=[1];
b2_num_A=[1]; b2_den_A=[L R];
b3_num_A=[-2*K*Io]; b3_den_A=[Zo^2];
b4_num_A=[1]; b4_den_A=[m];
b5_num_A=[1]; b5_den_A=[1 0 0];
%###########
%OLTF without compensator
[num1_A,den1_A]=series(b1_num_A,b1_den_A,b2_num_A,b2_den_A);
[num2_A,den2_A]=series(b3_num_A,b3_den_A,b4_num_A,b4_den_A);
[num3_A,den3_A]=series(num1_A,den1_A,num2_A,den2_A);
[GH_num_A, GH_den_A]=series(num3_A,den3_A,b5_num_A,b5_den_A);
unc_sys= tf(GH_num_A,GH_den_A);
%##################
%% PART A
m=100; % kg
ts=0.5; % 2%[s]
ksiw_n=4/(ts);
ksi=0.6;
des_pole=-ksiw_n+ksiw_n*tan(acos(ksi))*1i;
% Phase Difference and Desired Magnitude
phase_A=angle((polyval(GH_num_A,des_pole))/(polyval(GH_den_A,des_pole))); 
mag_A=abs((polyval(GH_num_A,des_pole))/(polyval(GH_den_A,des_pole)));
pole_angle_A=angle(des_pole); %beta
pole_mag_A=abs(des_pole);
rad2deg(pi-phase_A);
% Using Analytical PID Model for PD Controller
Kp_A=(-sin(pole_angle_A+phase_A))/(mag_A*sin(pole_angle_A));
Kd_A=(sin(phase_A))/(pole_mag_A*mag_A*sin(pole_angle_A));
% Compensator Block
Gc_num_A=[Kd_A Kp_A];
Gc_den_A=[1];
Gc_A=tf(Gc_num_A,Gc_den_A);
GH_C_num_A=conv(GH_num_A,Gc_num_A);
GH_C_den_A=conv(GH_den_A,Gc_den_A);
H=1; %Sensor TF
% OLTF
G_A=tf(GH_C_num_A,GH_C_den_A);
G_A=minreal(G_A);

% CLTF
CLTF_A=feedback(G_A,H);
CLTF_A=minreal(CLTF_A);
stepinfo(CLTF_A)
[CLTF_num_A,CLTF_den_A]=tfdata(CLTF_A);
% Transfer Functions of E/D and C/D with respectively for PD Controller 
s=tf('s');
ED_TF_A=(-(1/m)*(1/(s^2)))/(1+(-1*Gc_A*(1/(R+L*s))*((-2*K*Io)/(Zo^2))*(1/m)*(1/(s^2))));
ED_TF_A=minreal(ED_TF_A);
% Transfer Functions of C/D
CD_TF_A=((1/m)*(1/(s^2)))/(1+(-1*Gc_A*(1/(R+L*s))*((-2*K*Io)/(Zo^2))*(1/m)*(1/(s^2))));
CD_TF_A=minreal(CD_TF_A);
[CD_TF_num_A,CD_TF_den_A]=tfdata(CD_TF_A);
% Transfer Functions of E/R
ER_TF_A=1/(1+G_A)
%Response with Step Input and Step Disturbance
t=0:0.01:1.5;
% Step Input with Step Disturbances
stepA=step(CLTF_num_A{:},CLTF_den_A{:},t);
figure
plot(t,(stepA),'.-')
hold on
plot(t,ones(1,length(t)))
title('Step Response of PD Compensated Systems')
xlabel('Time[s]')
hold on
legend('Compensated System (A)','Input')
%##################
%% Part B(PID)
%System Properties
m=160; %total mass 
D1_B=m*9.81;
D2=-(K*Io^2)/(Zo^2);
%#############################
%Defining Blocks
b1_num_B=[inv]; b1_den_B=[1];
b2_num_B=[1]; b2_den_B=[L R];
b3_num_B=[-2*K*Io]; b3_den_B=[Zo^2];
b4_num_B=[1]; b4_den_B=[m];
b5_num_B=[1]; b5_den_B=[1 0 0];
%############################
%Finding OLTF
[num1_B,den1_B]=series(b1_num_B,b1_den_B,b2_num_B,b2_den_B);
[num2_B,den2_B]=series(b3_num_B,b3_den_B,b4_num_B,b4_den_B);
[num3_B,den3_B]=series(num1_B,den1_B,num2_B,den2_B);
[GH_num_B, GH_den_B]=series(num3_B,den3_B,b5_num_B,b5_den_B);
OLP_A=roots(GH_den_B);
%##################
phase_B=angle((polyval(GH_num_B,des_pole))/(polyval(GH_den_B,des_pole))); %phi
mag_B=abs((polyval(GH_num_B,des_pole))/(polyval(GH_den_B,des_pole)));

pole_angle_B=angle(des_pole); %beta
pole_mag_B=abs(des_pole);
Ki_B=(250*R*Zo^2*(D1_B+D2))/(K*Io);
Kp_B=((-sin(pole_angle_B+phase_B))/(mag_B*sin(pole_angle_B)))-((2*Ki_B*cos(pole_angle_B))/(pole_mag_B));
Kd_B=((sin(phase_B))/(pole_mag_B*mag_B*sin(pole_angle_B)))+(Ki_B/(pole_mag_B^2));
Gc_num_B=[Kd_B Kp_B Ki_B];
Gc_den_B=[1 0];
Gc_B=tf(Gc_num_B,Gc_den_B);

GH_c_num_B=conv(GH_num_B,Gc_num_B);
GH_c_den_B=conv(GH_den_B,Gc_den_B);
H=1; %Sensor TF
%OLTF
G_B=tf(GH_c_num_B,GH_c_den_B);
G_B=minreal(G_B);

CLTF_B=feedback(G_B,H);
CLTF_B=minreal(CLTF_B);
[CLTF_num_B,CLTF_den_B]=tfdata(CLTF_B);
roots(CLTF_den_B{:});
s=tf('s');
% E/D TF
ED_TF_B=(-s*(R+L*s)*Zo^2)/(((s*(R+L*s)*(Zo^2)*m*s^2))+(((Kd_B*s^2)+(Kp_B*s)+(Ki_B))*(2*K*Io)));
ED_TF_B=minreal(ED_TF_B)
% E/R TF
ER_TF_B=1/(1+G_B);
ER_TF_B=minreal(ER_TF_B)
% C/D TF for PID
CD_TF_B=((1/m)*(1/(s^2)))/(1+(-1*Gc_B*(1/(R+L*s))*((-2*K*Io)/(Zo^2))*(1/m)*(1/(s^2))));
CD_TF_B=minreal(CD_TF_B)
[CD_TF_num_B,CD_TF_den_B]=tfdata(CD_TF_B);

% Transient Response with Step Input and Step Disturbance
t=0:0.01:1.5;
stepB=step(CLTF_num_B{:},CLTF_den_B{:},t);
stepD_B=(D1_B+D2)*step(CD_TF_num_B{:},CD_TF_den_B{:},t);
figure
plot(t,(stepB+stepD_B),'.-')
hold on
plot(t,(1+stepD_B),'--')
hold on
axis([0 1.5 0 1.5]);
hold on
plot(t,ones(1,length(t)))
title('Step Response of PID Compensated System with Step Disturbances')
legend('Passenger Board Right at the Start','Passenger Board at Reference Value','Reference')
xlabel('Time[s]')
hold on
t=1.5:0.01:3.0;
stepB=step(CLTF_num_B{:},CLTF_den_B{:},t);
stepD_B=(D1_B+D2)*step(CD_TF_num_B{:},CD_TF_den_B{:},t);
plot(t,(stepB+stepD_B),'.-')
hold on
plot(t,(1+stepD_B),'--')
hold on
axis([0 1.5 0 1.75]);
grid on
hold on
plot(t,ones(1,length(t)))
title('Step Response of PID Compensated System with Step Disturbances')
legend('Passenger Board Right at the Start','Passenger Board at Reference Value','Reference')
xlabel('Time[s]')
hold off
%% PART C
%Response with Step Input and Step Disturbance
t1=0:0.01:1.5;
stepB_wo=step(CLTF_num_B{:},CLTF_den_B{:},t1);
t2=0:0.01:1.5;
stepD_B=(D1_B+D2)*step(CD_TF_num_B{:},CD_TF_den_B{:},t2);
figure
plot(t1,(stepB_wo),'.-')
hold on
last_val=0.9986;
plot(t,(last_val+stepD_B),'--');
grid on
axis([0 3 0 1.75]);
hold on
plot([t1+t2],ones(1,length(t1+t2)))
title('Step Response of PID Compensated System with Step Disturbances')
legend('Platform without the Passenger','Passenger Board at 1.5 Second','Reference')
xlabel('Time[s]')
hold off

% PART D
Required Info obtained in PART A and PART B
%% Part E case (a): without passenger.
clc;clear;close all
%#######################
% System Specifications
des_phase_margin=35; %degree
desired_gain_margin=14; %dB
%########################
%System Properties
N=100; %Turns
L=10; %H
R=500; %Ohm
Io=10; %A
m=100; %kg
A=20; %m^2
nu=2.9*10^-4; %F/m
g = 9.81; % m^2/s
K=(nu*A*(N^2))/2;  %Linearized Constant
Zo=1.72; %m
inv=-1; %AC->DC
%Defining Blocks
b1_num_a=[inv]; b1_den_a=[1];
b2_num_a=[1]; b2_den_a=[L R];
b3_num_a=[-2*K*Io]; b3_den_a=[Zo^2];
b4_num_a=[1]; b4_den_a=[m];
b5_num_a=[1]; b5_den_a=[1 0 0];
%############################
[num1_a,den1_a]=series(b1_num_a,b1_den_a,b2_num_a,b2_den_a);
[num2_a,den2_a]=series(b3_num_a,b3_den_a,b4_num_a,b4_den_a);
[num3_a,den3_a]=series(num1_a,den1_a,num2_a,den2_a);
[GH_num_a, GH_den_a]=series(num3_a,den3_a,b5_num_a,b5_den_a);
OLP_a=roots(GH_den_a);
%##################
%Lead Compensation
K=16000;
OLTF_uc=tf(GH_num_a,GH_den_a) % Uncompanstated OLTF
%##################
% bode(OLTF_uc*K);     %Bode Diagram
% nyquist(OLTF_uc*K);  %Nyquist Diagram
%margin(OLTF_uc*K);     %Margins of Uncompanstaed But Gain Asjusted Sytem
%##################
[mag_a,phase_a,wout_a]=bode(OLTF_uc*K);
mag_a=squeeze(mag_a);
phase_a=squeeze(phase_a);
[gain_margin_a,phase_margin_a,w_pc_a,w_gc_a]=margin(OLTF_uc*K);
phase_lead_a=des_phase_margin-phase_margin_a+8;
alpha_a=(1-sind(phase_lead_a))/(1+sind(phase_lead_a));
w_gc_new_a=interp1(mag_a,wout_a,sqrt(alpha_a));
T_a=1/(sqrt(alpha_a)*w_gc_new_a);
Kc_a=K/alpha_a;
Gc_num_a=[Kc_a*alpha_a*T_a Kc_a*alpha_a];
Gc_den_a=[alpha_a*T_a 1];
GH_compensated_num_a=conv(GH_num_a,Gc_num_a);
GH_compensated_den_a=conv(GH_den_a,Gc_den_a);
sys_compensated_a=tf(GH_compensated_num_a,GH_compensated_den_a);
[gain_margin_new_a,phase_margin_new_a,w_pc_new_A,w_gc_newnew_a]=margin(sys_compensated_a);
H=1; %Feedback TF
% Compansated OLTF
G_a=tf(GH_compensated_num_a,GH_compensated_den_a); 
G_a=minreal(G_a);

% CLTF
CLTF_a=feedback(G_a,H);
CLTF_a=minreal(CLTF_a)
[CLTF_num_a,CLTF_den_a]=tfdata(CLTF_a);
%Transient Behaviour 
t=0:0.01:1.5;
stepR=step(CLTF_num_a{:},CLTF_den_a{:},t);
figure
plot(t,(stepR),'.-')
grid on
hold on
plot(t,ones(1,length(t)))
title('Step Response of Lead Compensated System')
xlabel('Time[s]')
hold on
legend('Compensated System (a)','Input')
%% Part E case (a): with passenger.
clc;clear;close all;
%System Properties
inv=-1; %AC to DC
R=500; %Ohm
L=10; %H
nu=2.9*10^-4; %Farad/m
A=20; %m^2
N=100; %Turns
K=(nu*A*(N^2))/2;
Io=10; %A
Zo=1.72; %m
%#######################
% System Specifications
desired_phase_margin=35; %degree
desired_gain_margin=14; %dB
%#############################
m_b=160; %kg
D1_b=m_b*9.81;
D2=-(K*Io^2)/(Zo^2);
%Defining Blocks
b1_num_b=[inv]; b1_den_b=[1];
b2_num_b=[1]; b2_den_b=[L R];
b3_num_b=[-2*K*Io]; b3_den_b=[Zo^2];
b4_num_b=[1]; b4_den_b=[m_b];
b5_num_b=[1]; b5_den_b=[1 0 0];
%############################
%Finding OLTF
[num1_b,den1_b]=series(b1_num_b,b1_den_b,b2_num_b,b2_den_b);
[num2_b,den2_b]=series(b3_num_b,b3_den_b,b4_num_b,b4_den_b);
[num3_b,den3_b]=series(num1_b,den1_b,num2_b,den2_b);
[GH_num_b, GH_den_b]=series(num3_b,den3_b,b5_num_b,b5_den_b);
OLP_b=roots(GH_den_b);
%##################
%Lead Compensation
sys_b=tf(GH_num_b,GH_den_b);
K_b=34000;
sys_b=K_b*sys_b; %KG(s)
%##################
% bode(OLTF_uc*K);       %Bode Diagram
% nyquist(sys_b);        %Nyquist Diagram
% margin(OLTF_uc*K);     %Margins of Uncompanstaed But Gain Asjusted Sytem
%############
[mag_b,phase_b,wout_b]=bode(sys_b,{10^-4,10^4});
mag_b=squeeze(mag_b);
phase_b=squeeze(phase_b);
[gain_margin_b,phase_margin_b,w_pc_b,w_gc_b]=margin(sys_b); %Gain Adjusted but Uncompensated
phase_lead_needed_b=desired_phase_margin-phase_margin_b+12;
alpha_b=(1-sind(phase_lead_needed_b))/(1+sind(phase_lead_needed_b));
w_gc_found_b=interp1(mag_b,wout_b,sqrt(alpha_b));
T_b=1/(sqrt(alpha_b)*w_gc_found_b);
Kc_b=K_b/alpha_b;
Gc_num_b=[Kc_b*alpha_b*T_b Kc_b*alpha_b];
Gc_den_b=[alpha_b*T_b 1];
Gc_b=tf(Gc_num_b,Gc_den_b);
Gc_b=minreal(Gc_b);
GH_compensated_num_b=conv(GH_num_b,Gc_num_b);
GH_compensated_den_b=conv(GH_den_b,Gc_den_b);
G_b=tf(GH_compensated_num_b,GH_compensated_den_b);
G_b=minreal(G_b);
%Bode Plot of Gain Adjusted and Uncompensated System
figure
margin(G_b)
%############
[gain_margin_new_b,phase_margin_new_b,w_pc_new_b,w_gc_new_b]=margin(G_b); %Gain Adjusted and Compensated
H=1; %Sensor TF
CLTF_b=feedback(G_b,H);
CLTF_b=minreal(CLTF_b)
[CLTF_num_b,CLTF_den_b]=tfdata(CLTF_b);
CLP_b=roots(CLTF_den_b{:});
CLZ_b=roots(CLTF_num_b{:});
ER_TF_b=1/(1+G_b);
ER_TF_b=minreal(ER_TF_b);
s=tf('s');
ED_TF_b=(-(1/(m_b*s^2)))/(1+(Gc_b*((1/(R+L*s))*((2*K*Io)/(Zo^2)))*(1/(m_b*s^2))));
ED_TF_b=minreal(ED_TF_b)
CD_TF_b=((1/(m_b*s^2)))/(1+(Gc_b*((1/(R+L*s))*((2*K*Io)/(Zo^2)))*(1/(m_b*s^2))));
CD_TF_b=minreal(CD_TF_b);
[CD_TF_num_b,CD_TF_den_b]=tfdata(CD_TF_b);
%####################
%Responses
t_b=0:0.01:1.5;
t_a=0:0.01:1.5;
sa=step(CLTF_num_b{:},CLTF_den_b{:},t_a);
sb=step(CLTF_num_b{:},CLTF_den_b{:},t_b);
sD_b=(D1_b+D2)*step(CD_TF_num_b{:},CD_TF_den_b{:},t_b);
figure
plot(t_a,(sa),'-.')
hold on
plot([t_a 1.5+t_a],ones(1,2*length(t_a)),'-')
hold on
plot((1.5+t_b),(1+sD_b),'-.')
grid on
title('Unit Step Response Before and After the Mass Change')
xlabel('Time[s]')
legend('Platform without the Passenger','Reference','Passenger Board at 1.5 Second')
hold off







