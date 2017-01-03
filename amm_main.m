clear
global T T_ref T_pulse T_orig beta P V Q_in c_N2 c_H2 c_NH3 SDEN abyv c_tot pp
T_orig = 600;
T = T_orig;
T_pulse = 975;
P = 1;
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
Y_H2  = X_H2 /(X_H2+X_N2+X_NH3);
Y_N2  = X_N2 /(X_H2+X_N2+X_NH3);
Y_NH3 = X_NH3/(X_H2+X_N2+X_NH3);
V = 1;
Q_in = 3.76e-2;
SDEN = 3.1331e-10;             % Catalyst surface site density (moles/cm2)
abyv = 500;
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
options2 = odeset ('NonNegative',[1 2 3 4 5 6 7 8 9],'BDF','on','Stats','off','AbsTol',1e-12,'RelTol',1e-10);
tic;
%tspan = tspan1;
s0 = [0 0 0 0 0 0 c_N2 c_H2 c_NH3];
%[t,s] = ode15s(@ammonia,[0 tspan],s0,options0)
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
T_pulse = 600;
tspan = 1e6;%tspan1b - tspan;
[t,s] = ode15s(@ammonia,[0 tspan],s0,options2);
s0 = s(end,1:9);
%s(:,10) = SDEN - sum(s(:,2:6),2);
%save('grw3.mat','t','s','-v7.3');
clear t s;
tstart = 0;
pfrnodes = 1;
for pfr=1:pfrnodes;
T_pulse = 800;
tspan = 1000;%tspan2 - tspan;
[t,s] = ode15s(@ammonia,[tstart tspan+tstart],s0,options1);
s0 = zeros(1,9);
s0(7:9) = s(end,7:9);
s(:,10) = SDEN - sum(s(:,2:6),2);
tr{pfr}=t;
sr{pfr}=s;
tstart = t(end);
c_N2 = s0(7);
c_H2 = s0(8);
c_NH3 = s0(9);
end
%save('grw4.mat','t','s','-v7.3');
toc;
%save('ammonia_dual_pulse_2.mat','s','t')
figure(1)
hold on
for pfr=1:pfrnodes
plot(tr{pfr},sr{pfr}(:,7)./sum(sr{pfr}(:,7:9),2),'b')
plot(tr{pfr},sr{pfr}(:,8)./sum(sr{pfr}(:,7:9),2),'r')
plot(tr{pfr},sr{pfr}(:,9)./sum(sr{pfr}(:,7:9),2),'g')
end
hold off
legend('N_2','H_2','NH_3')
%xlim([tspan1 tspan2]);
fs=sr{pfrnodes}(end,7:9)./sum(sr{pfrnodes}(end,7:9),2)
figure(2)
hold on
for pfr=1:pfrnodes
plot(tr{pfr},sr{pfr}(:,1)./(SDEN),'b')
plot(tr{pfr},sr{pfr}(:,3)./(SDEN),'r')
plot(tr{pfr},sr{pfr}(:,4)./(SDEN),'g')
plot(tr{pfr},sr{pfr}(:,10)./(SDEN),'m')
end
hold off
legend('N_{2*}','H_*','NH_{3*}','\theta_*')
sr{pfrnodes}(end,[1:6 10])/(SDEN)