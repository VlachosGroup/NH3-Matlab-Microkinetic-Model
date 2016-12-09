clear
global T T_ref T_pulse T_orig beta P V SDEN abyv
T_orig = 675;
T = T_orig;
T_pulse = 975;
P = 200;
beta = 1;
tspan1 = 7472;
tspan1a = 1.5e4;
tspan2 = 5e5;
T_ref = 1; %K
X_H2  = 3;
X_N2  = 1;
X_NH3 = 0;
Y_H2 = X_H2/(X_H2+X_N2+X_NH3);
Y_N2 = X_N2/(X_H2+X_N2+X_NH3);
Y_NH3 = X_NH3/(X_H2+X_N2+X_NH3);
V = 1;
SDEN = 3.1331e-10;             % Catalyst surface site density (moles/cm2)
abyv = 800;
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
%options = odeset ('Jacobian',@dammonia,'NonNegative',[1 2 3 4 5 6 7 8 9],'Stats','on','AbsTol',1e-14,'RelTol',1e-12);
options1 = odeset ('MaxStep',0.001,'NonNegative',[1 2 3 4 5 6 7 8 9],'BDF','on','InitialStep',1e-10,'Stats','off','AbsTol',1e-14,'RelTol',1e-12);
options2 = odeset ('MaxStep',100,'NonNegative',[1 2 3 4 5 6 7 8 9],'BDF','on','InitialStep',1e-10,'Stats','off','AbsTol',1e-14,'RelTol',1e-12);
tic;
sol = ode15s(@ammonia,[0 tspan1],[0 0 0 0 0 0 c_N2 c_H2 c_NH3],options1);
s0 = sol.y';
T_pulse = 875;
sol  = odextend(sol,@ammonia,tspan1a,s0(end,:),options1);
s0 = sol.y';
sol  = odextend(sol,@ammonia,tspan2,s0(end,:),options2);
toc;
s = sol.y';
t = sol.x';
s(:,10) = SDEN - sum(s(:,2:6),2);
save('ammonia_dual_pulse_2.mat','s','t')
figure(1)
plot(t,s(:,[9 8 7])./sum(s(:,7:9),2),'-')
legend('NH_3','H_2','N_2')
%xlim([tspan1 tspan2]);
fs=s(end,7:9)./sum(s(end,7:9),2)
figure(2)
plot(t,s(:,[1 3 4 10])./(SDEN),'-')
legend('\theta_*','H_*','NH_*','N_{2*}')
s(end,[1:6 10])/(SDEN)