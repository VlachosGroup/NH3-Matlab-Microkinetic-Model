clear
global T T_ref T_pulse T_orig beta P V Q_in c_N2 c_H2 c_NH3 SDEN abyv c_tot pp
T_orig = 675;
T = T_orig;
T_pulse = 975;
P = 1.0;
pp=[0,P];
beta = 1;
tspan1 = 7473;
tspan1a = 14873;
tspan1b = 27000;
tspan2 = 5e5;
T_ref = 1; %K
X_H2  = 0;
X_N2  = 0;
X_NH3 = 1;
Y_H2 = X_H2/(X_H2+X_N2+X_NH3);
Y_N2 = X_N2/(X_H2+X_N2+X_NH3);
Y_NH3 = X_NH3/(X_H2+X_N2+X_NH3);
V = 1;
Q_in = 3.76e-3;
SDEN = 3.1331e-10;             % Catalyst surface site density (moles/cm2)
abyv = 100;
cat_surf = abyv*V;
MW_H = 1.007970;
MW_N = 14.006700;
R = 82.057338;  %(cm3 atm/K mol)
c_tot = P/(R*T);
c_H2 = Y_H2*c_tot;
c_N2 = Y_N2*c_tot;
c_NH3 = Y_NH3*c_tot;
M_H2 = c_H2*V*(2*MW_H);
M_N2 = c_N2*V*(2*MW_N);
M_NH3 = c_NH3*V*(3*MW_H+MW_N);
options0 = odeset ('MaxStep',0.001,'NonNegative',[1 2 3 4 5 6 7 8 9],'BDF','on','InitialStep',1e-10,'Stats','off','AbsTol',1e-14,'RelTol',1e-12);
options1 = odeset ('MaxStep',0.001,'NonNegative',[1 2 3 4 5 6 7 8 9],'BDF','on','InitialStep',1e-5,'Stats','off','AbsTol',1e-12,'RelTol',1e-10);
options2 = odeset ('MaxStep',0.1,'NonNegative',[1 2 3 4 5 6 7 8 9],'BDF','on','Stats','off','AbsTol',1e-12,'RelTol',1e-10);
tic;
%tspan = tspan1;
s0 = [0 0 0 0 0 0 c_N2 c_H2 c_NH3];
%[t,s] = ode15s(@ammonia,[0 tspan],s0,options0);
%s0 = s(end,1:9);
%s(:,10) = SDEN - sum(s(:,2:6),2);
%save('grw1.mat','t','s','-v7.3');
%clear t s;
T_pulse = 875;
%tspan = tspan1a - tspan;
%[t,s] = ode15s(@ammonia,[0 tspan],s0,options1);
%s0 = s(end,1:9);
%s(:,10) = SDEN - sum(s(:,2:6),2);
%save('grw2.mat','t','s','-v7.3');
%clear t s;
T_pulse = 675;
tspan = 1e3;%tspan1b - tspan;
[t,s] = ode15s(@ammonia,[0 tspan],s0,options2);
s0 = s(end,1:9);
%s(:,10) = SDEN - sum(s(:,2:6),2);
%save('grw3.mat','t','s','-v7.3');
clear t s;
T_pulse = 675;
tspan = 1e2;%tspan2 - tspan;
[t,s] = ode15s(@ammonia,[0 tspan],s0,options2);
s(:,10) = SDEN - sum(s(:,2:6),2);
%save('grw4.mat','t','s','-v7.3');
toc;
%save('ammonia_dual_pulse_2.mat','s','t')
figure(1)
plot(t,s(:,[9 8 7])./sum(s(:,7:9),2),'-')
legend('NH_3','H_2','N_2')
%xlim([tspan1 tspan2]);
fs=s(end,7:9)./sum(s(end,7:9),2)
figure(2)
plot(t,s(:,[1 3 4 10])./(SDEN),'-')
legend('\theta_*','H_*','NH_{3*}','N_{2*}')
s(end,[1:6 10])/(SDEN)