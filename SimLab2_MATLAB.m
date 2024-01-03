clear all; close all; clc
set(0, 'defaultLegendInterpreter','latex');
set(0, 'DefaultAxesTickLabelInterpreter','latex');
set(0,'defaultTextinterpreter','latex');

%% Part1 - 2nd Order Plant, 2nd Order Estimator
M = 20; m = 1; k = 32; b = 0.3; f = 5; T = 1/f; 
y0 = 0; ydot0 = 0; d0 = 0; ddot0 = 0; xhat0 = 0;
sat = 10; noise = 0; Tf = 20;
F = [0 0; 1 0]; G = [1/(M+m);0]; H = [0 1]; J = 0;
sysc = ss(F,G,H,J);
[sysd,~] = c2d(sysc,T,'zoh');
[Phi,Gamma,H,J] = ssdata(sysd)
G2 = tf(sysd)

%%
% 1a
C = [Gamma Phi*Gamma]; %ctrb(Phi,Gamma)
detC = det(C)
wn = 1; zeta = 0.5;
s1 = -zeta*wn + j*wn*sqrt(1-zeta^2);
s2 = -zeta*wn - j*wn*sqrt(1-zeta^2);
pole1 = exp(s1*T); pole2 = exp(s2*T);
K = place(Phi,Gamma,[pole1, pole2])
NxNu = [Phi-eye(2) Gamma; H J]\[zeros(2,1);1]; % inv(A)*b
Nx = NxNu(1:end-1); Nu = NxNu(end);
N = Nu + K*Nx

%%
% 1b
O = [H; H*Phi]; %obsv(Phi,H)
detO = det(O)
wn = 4; zeta = 0.7;
s1 = -zeta*wn + j*wn*sqrt(1-zeta^2);
s2 = -zeta*wn - j*wn*sqrt(1-zeta^2);
p1 = exp(s1*T); p2 = exp(s2*T);
L = place(Phi',H',[p1, p2])'
D2 = ss(Phi-Gamma*K-L*H, L, -K, [], T);
D2 = tf(D2)
h22 = feedback(G2,-D2); H22 = N*h22;

%%
set_param([gcs,'/ManualSwitch5'],'sw','0')
set_param([gcs,'/ReducedEstimator'],'Commented','on')
set_param([gcs,'/K1'],'Commented','on')
set_param([gcs,'/Estimator'],'Commented','off')
set_param([gcs,'/K'],'Commented','off')
% set_param([gcs,'/Subtract'],'Commented','off')
% set_param([gcs,'/Scope'],'Commented','off')
% set_param([gcs,'/To Workspace1'],'Commented','off')
set_param([gcs,'/Mux'],'Commented','on')
set_param([gcs,'/filter'],'Commented','off')
set_param([gcs,'/ManualSwitch'],'sw','0')
set_param([gcs,'/ManualSwitch1'],'sw','0')
set_param([gcs,'/ManualSwitch2'],'sw','1')
set_param([gcs,'/ManualSwitch3'],'sw','1')
set_param([gcs,'/noise/ManualSwitch4'],'sw','0')
set_param([gcs,'/noise1/ManualSwitch4'],'sw','0')
set_param([gcs,'/2nd-Order Plant'],'Commented','off')
set_param([gcs,'/4th-Order Plant'],'Commented','on')
set_param([gcs,'/Saturation'],'Commented','through') % without
set_param([gcs,'/Saturation1'],'Commented','through')
out = sim('SimulationLab2.slx',Tf);
set_param([gcs,'/Saturation'],'Commented','off') % with
set_param([gcs,'/Saturation1'],'Commented','off')
out2 = sim('SimulationLab2.slx',Tf);

%%
%a
figure();
plot(out.r,'r',"LineWidth",1.5); hold on
u_max = round(max(out.u));
plot(out.u/u_max,'-.b',"LineWidth",1.5);%stairs(out.time.Data,out.u.Data/u_max,'b');
hold on;
u_max2 = round(max(out2.u)); 
plot(out2.u/u_max2,'--g',"LineWidth",1.5); hold on
plot(out.y,'-.k',"LineWidth",1.5); hold on
plot(out2.y,'--m',"LineWidth",1.5);
legend("r","u/"+num2str(u_max)+" -without saturation",...
            "u/"+num2str(u_max2)+" -with saturation",...
            "y -without saturation","y -with saturation",...
            'FontSize',12,"Location","best")
xlabel("t(s)"); ylabel(""); grid on
title("\textbf{1a) Plots for full state feedback}",'interpreter','latex')

%b
figure();
plot(out.xtilde,"LineWidth",1.5); %'-.b',
hold on; 
plot(out2.xtilde,"LineWidth",1.5); hold on %'--g',
legend("$\tilde{\dot{y}}$ -without saturation",...
         "$\tilde{y}$ -without saturation",...
         "$\tilde{\dot{y}}$ -with saturation",...
         "$\tilde{y}$ -with saturation",...
         'FontSize',12,"Location","best")
xlabel("t(s)"); ylabel(""); grid on
title("\textbf{1b) Plots of $\bf \tilde{x}$}",'interpreter','latex')

%%
% 1c
%set_param([gcs,'/Subtract'],'Commented','on')
%set_param([gcs,'/Scope'],'Commented','on')
%set_param([gcs,'/To Workspace1'],'Commented','on')
set_param([gcs,'/ManualSwitch2'],'sw','0')
set_param([gcs,'/Saturation'],'Commented','through') % without
set_param([gcs,'/Saturation1'],'Commented','through')
out3 = sim('SimulationLab2.slx',Tf);
set_param([gcs,'/Saturation'],'Commented','off') % with
set_param([gcs,'/Saturation1'],'Commented','off')
out4 = sim('SimulationLab2.slx',Tf);

%%
figure();
plot(out3.r,'r',"LineWidth",1.5); hold on
u_max = round(max(out3.u));
plot(out3.u/u_max,'-.b',"LineWidth",1.5);
hold on;
u_max2 = round(max(out4.u)); 
plot(out4.u/u_max2,'--g',"LineWidth",1.5); hold on
plot(out3.y,'-.k',"LineWidth",1.5); hold on
plot(out4.y,'--m',"LineWidth",1.5);
legend("r","u/"+num2str(u_max)+" -without saturation",...
            "u/"+num2str(u_max2)+" -with saturation",...
            "y -without saturation","y -with saturation",...
            'FontSize',12,"Location","best")
xlabel("t(s)"); ylabel(""); grid on
title("\textbf{1c) Plots for combined control law and estimator}",'interpreter','latex')

%% 
% 1d
%with step disturbance
noise = 0.1;
set_param([gcs,'/noise/ManualSwitch4'],'sw','0')
set_param([gcs,'/noise1/ManualSwitch4'],'sw','0')
out5 = sim('SimulationLab2.slx',Tf);

%with random disturbance
noise = 0.01;
set_param([gcs,'/noise/ManualSwitch4'],'sw','1')
set_param([gcs,'/noise1/ManualSwitch4'],'sw','1')
out6 = sim('SimulationLab2.slx',Tf);

%%
figure();
plot(out5.r,'r',"LineWidth",1.5); hold on
u_max = round(max(out5.u));
plot(out5.u/u_max,'-.b',"LineWidth",1.5);
hold on;
u_max2 = round(max(out6.u)); 
plot(out6.u/u_max2,'--g',"LineWidth",1.5); hold on
plot(out5.y,'-.k',"LineWidth",1.5); hold on
plot(out6.y,'--m',"LineWidth",1.5);
legend("r","u/"+num2str(u_max)+" -with step disturbance",...
            "u/"+num2str(u_max2)+" -with random disturbance",...
            "y -with step disturbance","y -with random disturbance",...
            'FontSize',12,"Location","best")
xlabel("t(s)"); ylabel(""); grid on
title("\textbf{1d) Plots for combined compensator with noise}",'interpreter','latex')

%%
% <html>
% <p style="page-break-before: always">
% </html>

%% Part2 - 4th Order Plant, 2nd Order Estimator
%Without Noise
noise = 0;
set_param([gcs,'/2nd-Order Plant'],'Commented','on')
set_param([gcs,'/4th-Order Plant'],'Commented','off')
set_param([gcs,'/ManualSwitch'],'sw','1')
set_param([gcs,'/ManualSwitch1'],'sw','1')
set_param([gcs,'/ManualSwitch2'],'sw','0')
set_param([gcs,'/noise/ManualSwitch4'],'sw','0')
set_param([gcs,'/noise1/ManualSwitch4'],'sw','0')
set_param([gcs,'/Saturation'],'Commented','through') % without
set_param([gcs,'/Saturation1'],'Commented','through')
out7 = sim('SimulationLab2.slx',Tf);
set_param([gcs,'/Saturation'],'Commented','off') % with
set_param([gcs,'/Saturation1'],'Commented','off')
out8 = sim('SimulationLab2.slx',Tf);

%%
figure();
plot(out7.r,'r',"LineWidth",1.5); hold on
u_max = round(max(out7.u));
plot(out7.u/u_max,'-.b',"LineWidth",1.5);
hold on;
u_max2 = round(max(out8.u)); 
plot(out8.u/u_max2,'--g',"LineWidth",1.5); hold on
plot(out7.y,'-.k',"LineWidth",1.5); hold on
plot(out8.y,'--m',"LineWidth",1.5);
legend("r","u/"+num2str(u_max)+" -without saturation",...
            "u/"+num2str(u_max2)+" -with saturation",...
            "y -without saturation","y -with saturation",...
            'FontSize',12,"Location","best")
xlabel("t(s)"); ylabel(""); grid on
title("\textbf{2) Plots for combined control law and estimator}",'interpreter','latex')

%%
%With noise and saturation (2,3e)
%with step disturbance
noise = 0.1;
set_param([gcs,'/noise/ManualSwitch4'],'sw','0')
set_param([gcs,'/noise1/ManualSwitch4'],'sw','0')
out9 = sim('SimulationLab2.slx',Tf);
% with random disturbance
noise = 0.01;
set_param([gcs,'/noise/ManualSwitch4'],'sw','1')
set_param([gcs,'/noise1/ManualSwitch4'],'sw','1')
out10 = sim('SimulationLab2.slx',Tf);

%%
figure();
plot(out9.r,'r',"LineWidth",1.5); hold on
u_max = round(max(out9.u));
plot(out9.u/u_max,'-.b',"LineWidth",1.5);
hold on;
u_max2 = round(max(out10.u)); 
plot(out10.u/u_max2,'--g',"LineWidth",1.5); hold on
plot(out9.y,'-.k',"LineWidth",1.5); hold on
plot(out10.y,'--m',"LineWidth",1.5);
legend("r","u/"+num2str(u_max)+" -with step disturbance",...
            "u/"+num2str(u_max2)+" -with random disturbance",...
            "y -with step disturbance","y -with random disturbance",...
            'FontSize',8,"Location","best")
xlabel("t(s)"); ylabel(""); grid on
title("\textbf{2) Plots for combined compensator with noise}",'interpreter','latex')

%%
% <html>
% <p style="page-break-before: always">
% </html>

%% Part3 - 4th Order Plant, 4th Order Estimator
F = [-b/M -k/M b/M k/M; 
     1 0 0 0;
     b/m k/m -b/m -k/m;
     0 0 1 0];
G = [1/M; 0; 0; 0];
H = [0 0 0 1];
J = 0;
sysc = ss(F,G,H,J);
[sysd,~] = c2d(sysc,T,'zoh');
[Phi,Gamma,H,J] = ssdata(sysd)
G4 = tf(sysd)
h42 = feedback(G4,-D2); H42 = N*h42;

%%
% 3a
noise = 0;
zeta = 0.3; wn = sqrt((1/m + 1/M)*k);
P = zeros(1,4);
P(1) = pole1; P(2) = pole2; % from 1a
s3 = -zeta*wn + j*wn*sqrt(1-zeta^2);
s4 = -zeta*wn - j*wn*sqrt(1-zeta^2);
P(3) = exp(s3*T);
P(4) = exp(s4*T)
C = [Gamma Phi*Gamma Phi^2*Gamma Phi^3*Gamma];
detC = det(C) % not 0
K = place(Phi,Gamma,P)
NxNu = [Phi-eye(4) Gamma; H J]\[zeros(4,1);1]; % inv(A)*b
Nx = NxNu(1:end-1); Nu = NxNu(end);
N = Nu + K*Nx

%%
% 3b
O = [H; H*Phi; H*Phi^2; H*Phi^3];
detO = det(O)
Po = P/2.5;
L = place(Phi',H',Po)'
D4 = ss(Phi-Gamma*K-L*H, L, -K, [], T);
D4 = tf(D4)
h44 = feedback(G4,-D4); H44 = N*h44;

%%
set_param([gcs,'/Mux'],'Commented','off')
set_param([gcs,'/filter'],'Commented','on')
set_param([gcs,'/ManualSwitch2'],'sw','1')
set_param([gcs,'/ManualSwitch3'],'sw','0')
set_param([gcs,'/noise/ManualSwitch4'],'sw','0')
set_param([gcs,'/noise1/ManualSwitch4'],'sw','0')
set_param([gcs,'/2nd-Order Plant'],'Commented','on')
set_param([gcs,'/4th-Order Plant'],'Commented','off')
set_param([gcs,'/Saturation'],'Commented','through') % without
set_param([gcs,'/Saturation1'],'Commented','through')
out11 = sim('SimulationLab2.slx',Tf);
set_param([gcs,'/Saturation'],'Commented','off') % with
set_param([gcs,'/Saturation1'],'Commented','off')
out12 = sim('SimulationLab2.slx',Tf);

%%
%a
figure();
plot(out11.r,'r',"LineWidth",1.5); hold on
u_max = round(max(out11.u));
plot(out11.u/u_max,'-.b',"LineWidth",1.5);
hold on;
u_max2 = round(max(out12.u)); 
plot(out12.u/u_max2,'--g',"LineWidth",1.5); hold on
plot(out11.y,'-.k',"LineWidth",1.5); hold on
plot(out12.y,'--m',"LineWidth",1.5);
legend("r","u/"+num2str(u_max)+" -without saturation",...
            "u/"+num2str(u_max2)+" -with saturation",...
            "y -without saturation","y -with saturation",...
            'FontSize',12,"Location","best")
xlabel("t(s)"); ylabel(""); grid on
title("\textbf{3a) Plots for full state feedback}",'interpreter','latex')

%b
figure();
plot(out11.xtilde,"LineWidth",1.5); %'-.b',
hold on; 
plot(out12.xtilde,"LineWidth",1.5); hold on %'--g',
legend("$\tilde{\dot{y}}$ -without saturation",...
         "$\tilde{y}$ -without saturation",...
         "$\tilde{\dot{d}}$ -without saturation",...
         "$\tilde{d}$ -without saturation",...
          "$\tilde{\dot{y}}$ -with saturation",...
          "$\tilde{y}$ -with saturation",...
          "$\tilde{\dot{d}}$ -with saturation",...
          "$\tilde{d}$ -with saturation",...
          'FontSize',10,"Location","best")
xlabel("t(s)"); ylabel(""); grid on
title("\textbf{3b) Plots of $\bf \tilde{x}$}",'interpreter','latex')

%%
% 3c
set_param([gcs,'/ManualSwitch2'],'sw','0')
set_param([gcs,'/Saturation'],'Commented','through') % without
set_param([gcs,'/Saturation1'],'Commented','through')
out13 = sim('SimulationLab2.slx',Tf);
set_param([gcs,'/Saturation'],'Commented','off') % with
set_param([gcs,'/Saturation1'],'Commented','off')
out14 = sim('SimulationLab2.slx',Tf);

%%
figure();
plot(out13.r,'r',"LineWidth",1.5); hold on
u_max = round(max(out13.u));
plot(out13.u/u_max,'-.b',"LineWidth",1.5);
hold on;
u_max2 = round(max(out14.u)); 
plot(out14.u/u_max2,'--g',"LineWidth",1.5); hold on
plot(out13.y,'-.k',"LineWidth",1.5); hold on
plot(out14.y,'--m',"LineWidth",1.5);
legend("r","u/"+num2str(u_max)+" -without saturation",...
            "u/"+num2str(u_max2)+" -with saturation",...
            "y -without saturation","y -with saturation",...
            'FontSize',12,"Location","best")
xlabel("t(s)"); ylabel(""); grid on
title("\textbf{3c) Plots for combined control law and estimator}",'interpreter','latex')

%%
% 3d
init = 2;
ydot0 = init; y0 = 0; ddot0 = 0; d0 = 0;
out15 = sim('SimulationLab2.slx',Tf+5);
ydot0 = 0; y0 = init; ddot0 = 0;  d0 = 0; 
out16 = sim('SimulationLab2.slx',Tf+5);
ydot0 = 0; y0 = 0; ddot0 = init; d0 = 0;
out17 = sim('SimulationLab2.slx',Tf+5);
ydot0 = 0; y0 = 0; ddot0 = 0; d0 = init;
out18 = sim('SimulationLab2.slx',Tf);
d0 = 0; xhat0 = [init,init+1,init-1,init]; %ydot,y,ddot,d
out19 = sim('SimulationLab2.slx',Tf);
ydot0 = round(rand(1)*100)/100; y0 = round(rand(1)*100)/100;
ddot0 = round(rand(1)*100)/100; d0 = round(rand(1)*100)/100;
xhat0 = round([2*rand(1),2*rand(1),2*rand(1),2*rand(1)]*100)/100;
out20 = sim('SimulationLab2.slx',Tf+5);
x0 = ["2,0,0,0", "0,2,0,0", "0,0,2,0", "0,0,0,2","0,0,0,0",...
        regexprep(num2str([ydot0,y0,ddot0,d0]),'\s+',',')];
xHat0 = ["0,0,0,0","0,0,0,0","0,0,0,0","0,0,0,0","2,3,1,2",...
            regexprep(num2str(xhat0),'\s+',',')];
out = [out15,out16,out17,out18,out19,out20];

for i = 1:6
    figure();
    plot(out(i).r,'r',"LineWidth",1.5); hold on
    u_max = round(max(out(i).u));
    plot(out(i).u/u_max,'-.b',"LineWidth",1.5); hold on;
    plot(out(i).y,'-.k',"LineWidth",1.5);
    legend("r","u/"+num2str(u_max),"y",...
                'FontSize',10,"Location","best")
    xlabel("t(s)"); ylabel(""); grid on
    title("3d) Initial Conditions: x = ["+x0(i)+...
        "], $\hat{x}$= ["+xHat0(i)+"]",'interpreter','latex')  
end
y0 = 0; ydot0 = 0; d0 = 0; ddot0 = 0; xhat0 = 0;

%%
% 3e
%with noise and saturation
%with step disturbance
noise = 0.1;
set_param([gcs,'/noise/ManualSwitch4'],'sw','0')
set_param([gcs,'/noise1/ManualSwitch4'],'sw','0')
out21 = sim('SimulationLab2.slx',Tf);
% with random disturbancerank(ans
noise = 0.01;
set_param([gcs,'/noise/ManualSwitch4'],'sw','1')
set_param([gcs,'/noise1/ManualSwitch4'],'sw','1')
out22 = sim('SimulationLab2.slx',Tf);

%%
%compensator
figure()
rlocus(-D2*G4); hold on
axis([-8.2,2,-1.5,1.5])
r42 = rlocus(-D2*G4,1);
scatter(real(r42),imag(r42),50,"k","s","filled");
title("3e) Root locus with 2nd Order compensator")
legend("Root Locus","Closed Loop Poles",'Location','best')
set(gca,'XGrid','on')
set(gca,'YGrid','on')

figure()
rlocus(-D4*G4); hold on;
axis([-8.2,2,-1.5,1.5])
r44 = rlocus(-D4*G4,1);
scatter(real(r44),imag(r44),50,"k","s","filled");
title("3e) Root locus with 4th Order compensator")
legend("Root Locus","Closed Loop Poles",'Location','best')
set(gca,'XGrid','on')
set(gca,'YGrid','on')

% noise stuff
figure();
plot(out21.r,'r',"LineWidth",1.5); hold on
u_max = round(max(out21.u));
plot(out21.u/u_max,'-.b',"LineWidth",1.5);
hold on;
u_max2 = round(max(out22.u)); 
plot(out22.u/u_max2,'--g',"LineWidth",1.5); hold on
plot(out21.y,'-.k',"LineWidth",1.5); hold on
plot(out22.y,'--m',"LineWidth",1.5);
legend("r","u/"+num2str(u_max)+" -with step disturbance",...
            "u/"+num2str(u_max2)+" -with random disturbance",...
            "y -with step disturbance","y -with random disturbance",...
            'FontSize',8,"Location","best")
xlabel("t(s)"); ylabel(""); grid on
title("\textbf{3e) Plots for combined compensator with noise}",'interpreter','latex')

%%
% 3f
noise = 0;
% changed (lowered) f but was sensitive
zeta = 0.8; wn = 0.6*sqrt((1/m + 1/M)*k);
P = zeros(1,4);
P(1) = 0.98*pole1; P(2) = 0.98*pole2; %updated from zeta = 0.5, wn = 1
s3 = -zeta*wn + j*wn*sqrt(1-zeta^2);
s4 = -zeta*wn - j*wn*sqrt(1-zeta^2);
P(3) = exp(s3*T);
P(4) = exp(s4*T)
C = [Gamma Phi*Gamma Phi^2*Gamma Phi^3*Gamma];
detC = det(C); % not 0
K = place(Phi,Gamma,P)
NxNu = [Phi-eye(4) Gamma; H J]\[zeros(4,1);1]; % inv(A)*b
Nx = NxNu(1:end-1); Nu = NxNu(end);
N = Nu + K*Nx

O = [H; H*Phi; H*Phi^2; H*Phi^3];
detO = det(O);
Po = P/2.5;
L = place(Phi',H',Po)'
D4 = ss(Phi-Gamma*K-L*H, L, -K, [], T);
D4 = tf(D4)

set_param([gcs,'/Saturation'],'Commented','through') % without
set_param([gcs,'/Saturation1'],'Commented','through')
out25 = sim('SimulationLab2.slx',Tf);
set_param([gcs,'/Saturation'],'Commented','off') % with
set_param([gcs,'/Saturation1'],'Commented','off')
out26 = sim('SimulationLab2.slx',Tf);
figure();
plot(out25.r,'r',"LineWidth",1.5); hold on
u_max = round(max(out25.u));
plot(out25.u/u_max,'-.b',"LineWidth",1.5); hold on;
u_max2 = round(max(out26.u)); 
plot(out26.u/u_max2,'--g',"LineWidth",1.5); hold on
plot(out25.y,'-.k',"LineWidth",1.5); hold on
plot(out26.y,'--m',"LineWidth",1.5);
legend("r","u/"+num2str(u_max)+" -without saturation",...
            "u/"+num2str(u_max2)+" -with saturation",...
            "y -without saturation","y -with saturation",...
            'FontSize',8,"Location","best")
xlabel("t(s)"); ylabel(""); grid on
title("\textbf{3f) Unsaturated Example}",'interpreter','latex')

zeta = 0.3; wn = sqrt((1/m + 1/M)*k);
P = zeros(1,4);
P(1) = pole1; P(2) = pole2; % from 1a
s3 = -zeta*wn + j*wn*sqrt(1-zeta^2);
s4 = -zeta*wn - j*wn*sqrt(1-zeta^2);
P(3) = exp(s3*T);
P(4) = exp(s4*T);
C = [Gamma Phi*Gamma Phi^2*Gamma Phi^3*Gamma];
detC = det(C);
K = place(Phi,Gamma,P);
NxNu = [Phi-eye(4) Gamma; H J]\[zeros(4,1);1]; % inv(A)*b
Nx = NxNu(1:end-1); Nu = NxNu(end);
N = Nu + K*Nx;
% 3b
O = [H; H*Phi; H*Phi^2; H*Phi^3];
detO = det(O);
Po = P/2.5;
L = place(Phi',H',Po)';
D4 = ss(Phi-Gamma*K-L*H, L, -K, [], T);
D4 = tf(D4);

%%
% <html>
% <p style="page-break-before: always">
% </html>

%% Part4 - Bonus
%a - too much code handling
%%
%b - too many parameter options
%%
%c
noise = 0;
Phi11 = Phi(end,end); Phi12 = Phi(end,1:end-1); % d is last row
Phi21 = Phi(1:end-1,end); Phi22 = Phi(1:end-1,1:end-1);
Gamma1 = Gamma(end); Gamma2 = Gamma(1:end-1);
Pr = [0,1e-4,-1e-4];
Lr = place(Phi22',Phi12',Pr)'
%Phi21LrPhi11 = Phi21 - Lr*Phi11;

set_param([gcs,'/ManualSwitch5'],'sw','1')
set_param([gcs,'/ReducedEstimator'],'Commented','off')
set_param([gcs,'/K1'],'Commented','off')
set_param([gcs,'/Estimator'],'Commented','on')
set_param([gcs,'/K'],'Commented','on')
set_param([gcs,'/Saturation'],'Commented','through') % without
set_param([gcs,'/Saturation1'],'Commented','through')
out23 = sim('SimulationLab2.slx',Tf);
set_param([gcs,'/Saturation'],'Commented','off') % with
set_param([gcs,'/Saturation1'],'Commented','off')
out24 = sim('SimulationLab2.slx',Tf);

figure();
plot(out23.r,'r',"LineWidth",1.5); hold on
u_max = round(max(out23.u));
plot(out23.u/u_max,'-.b',"LineWidth",1.5);
hold on;
u_max2 = round(max(out24.u)); 
plot(out24.u/u_max2,'--g',"LineWidth",1.5); hold on
plot(out23.y,'-.k',"LineWidth",1.5); hold on
plot(out24.y,'--m',"LineWidth",1.5);
legend("r","u/"+num2str(u_max)+" -without saturation",...
            "u/"+num2str(u_max2)+" -with saturation",...
            "y -without saturation","y -with saturation",...
            'FontSize',12,"Location","best")
xlabel("t(s)"); ylabel(""); grid on
title("\textbf{4c) Plots with reduced state estimator}",'interpreter','latex')

%%
%d
opts = bodeoptions;
opts.MagUnits = 'abs'; opts.MagScale = 'log';
figure()
bodeplot(-D2,opts); grid on
title("4d) Bode Diagram of 2nd Order Compensator D_2")

figure()
bodeplot(-D4,opts); grid on
title("4d) Bode Diagram of 4th Order Compensator D_4")

figure()
bodeplot(-D2*G2,opts); grid on
title("4d) Bode Diagram of D_2G_2")

figure()
bodeplot(-D2*G4,opts); grid on
title("4d) Bode Diagram of D_2G_4")

figure()
bodeplot(-D4*G4,opts); grid on
title("4d) Bode Diagram of D_4G_4")

figure()
bodeplot(H22,opts); grid on
title("4d) Bode Diagram of H_{22}")

figure()
bodeplot(H42,opts); grid on
title("4d) Bode Diagram of H_{42}")

figure()
bodeplot(H44,opts); grid on
title("4d) Bode Diagram of H_{44}")
