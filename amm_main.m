%%       -------------------------------------------------
%                     NH3 Micro-kinetic model
%
%          Gerhard R Wittreich, P.E.  (February 10, 2017)
%        --------------------------------------------------
%
%  Main program: 
%      requires: ammonia.m, amm_kinetics.m, amm_thermo.m, amm_coverage.m
%                amm_BEP_LSR.m
%%
clear
datetime
%  Set key model parameters
%
global T T_ref T_pulse T_orig beta P V Q_in c_N2 c_H2 c_NH3 SDEN abyv ...
       c_tot Stoic_surf Stoic_gas Stoic MWON Isobaric Ea A Stick R_e ...
       R_k R MW_N2 MW_H2 MW_NH3
T_orig = 700;             % Reactor bulk temperature [K]
T = T_orig;               % Initial reactor temperature [K]
P = 1;                    % Reactor pressure [atm]
beta = [0 1 0 1 1 1 0]';
Ea   = [0 0 0 0 0 0 0]';
A    = [1.16e18 2.05e19 1.06e20 8.38e19]';
Stick= [0.5 0.5 0.5]';
MWON = 0;           % Motz-Wise Correction: 0 = Off (S); 1 = On (S/(1-S/2))
Isobaric = 1;       % Isobaric Reactor (1=Isobaric, 2=Isochoric)
T_ref = 1;          % Reference Temperature [K]
MW_H = 1.007970;                %---
MW_N = 14.006700;               %
MW_H2 = 2*MW_H;                 % Molecular mass [g/mol]
MW_N2 = 2*MW_N;                 %
MW_NH3 = MW_N + 3*MW_H;         %---
X_H2  = 0;                      %
X_N2  = 0;                      % Initial mole fraction in reactor feed
X_NH3 = 1;                      %
Y_H2  = X_H2 /(X_H2+X_N2+X_NH3);% Mole fractions
Y_N2  = X_N2 /(X_H2+X_N2+X_NH3);% normalized
Y_NH3 = X_NH3/(X_H2+X_N2+X_NH3);% to 1
V = 10;                         % Reactor volume (cm3)
Q_in = 3.76e-2;                % 0 = Batch Reactor,  Any other value = CSTR
SDEN = 4.4385e-10;             % Catalyst surface site density (moles/cm2)
abyv = 3000;                   % Catalyst loading (cm2 catalyst/cm3 reac volume)
R_e = 1.987e-3;                % Gas constant, (kcal/mol K)
R_k = 8.31451e7;               % Gas constant, (g cm2/mol K s)
R = 82.057338;                 % Gas constant, (cm3 atm/K mol)
c_tot = P/(R*T);               % Total starting moles in reactor [mol]
c_H2 = Y_H2*c_tot;             % Moles H2 in reactor [mol]
c_N2 = Y_N2*c_tot;             % Moles N2 in reactor [mol]
c_NH3 = Y_NH3*c_tot;           % Moles NH3 in reactor [mol]
M_H2 = c_H2*V*MW_H2;           % Mass H2 in reactor [kg]
M_N2 = c_N2*V*MW_N2;           % Mass N2 in reactor [kg]
M_NH3 = c_NH3*V*MW_NH3;        % Mass NH3 in reactor [kg]
Stoic_surf = [ 1  0  0  0  0  0  0  0  0 -1;... % Reaction
              -1  2  0  0  0  0  0  0  0 -1;... %
               0  0  2  0  0  0  0  0  0 -2;... % Surface
               0  0  1 -1  1  0  0  0  0 -1;... %
               0  0  1  0 -1  1  0  0  0 -1;... % Stoichiometry
               0  1  1  0  0 -1  0  0  0 -1;... %
               0  0  0  1  0  0  0  0  0 -1];   % ---
Stoic_gas =  [ 0  0  0  0  0  0 -1  0  0  0;... % Reaction
               0  0  0  0  0  0  0  0  0  0;... %
               0  0  0  0  0  0  0 -1  0  0;... % Gas
               0  0  0  0  0  0  0  0  0  0;... %
               0  0  0  0  0  0  0  0  0  0;... % Stoichiometry
               0  0  0  0  0  0  0  0  0  0;... %
               0  0  0  0  0  0  0  0 -1  0];   % ---
Stoic = Stoic_surf + Stoic_gas;                 % Total stoichiometry
% ODE Solver options
options0 = odeset ('MaxStep',0.001,'NonNegative',[1 2 3 4 5 6 7 8 9],...
                   'BDF','off','InitialStep',1e-10,'Stats','off',...
                   'AbsTol',1e-14,'RelTol',1e-12);
options1 = odeset ('MaxStep',0.0001,'NonNegative',[1 2 3 4 5 6 7 8 9],...
                   'BDF','off','InitialStep',1e-10,'Stats','off',...
                   'AbsTol',1e-14,'RelTol',1e-12);
options2 = odeset ('NonNegative',[1 2 3 4 5 6 7 8 9],'InitialStep',1e-10,...
                   'BDF','on','Stats','off','AbsTol',1e-14,'RelTol',1e-12);
tic;
s0 = [0 0 0 0 0 0 c_N2 c_H2 c_NH3 SDEN*abyv];             % Initial species concentrations
if ne(1,1)
    T_pulse = T_orig;
    tspan = max(floor(10/Q_in),5);
    [t,s] = ode15s(@ammonia,[0 tspan],s0,options2);
    s0 = s(end,1:9);
    clear t s
    %s0(:,10) = (SDEN*abyv) - sum(s(:,1:6),2);
end
tstart = 0;
pfrnodes = 1;           % PFR capability is not implemented.  Must be 1.
for pfr=1:pfrnodes;
    T_pulse = 673;
    tspan = max(floor(3/Q_in)+.5,2);
    [t,s] = ode15s(@ammonia,[tstart tspan+tstart],s0,options2);
    s(:,10) = (SDEN*abyv) - sum(s(:,1:6),2);
    tr{pfr}=t;
    sr{pfr}=s;
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
fs = sr{pfrnodes}(end,7:9)./sum(sr{pfrnodes}(end,7:9),2);
fprintf('\n------------------------\n')
fprintf('      Gas Species\n')
fprintf('------------------------\n')
fprintf('  N2       H2       NH3\n')
fprintf('------   ------   ------\n')
fprintf('%6.4f   %6.4f   %6.4f\n',fs)
figure(2)
hold on
for pfr=1:pfrnodes
plot(tr{pfr},sr{pfr}(:,1)./(SDEN*abyv),'b')
plot(tr{pfr},sr{pfr}(:,2)./(SDEN*abyv),'r')
plot(tr{pfr},sr{pfr}(:,3)./(SDEN*abyv),'c')
plot(tr{pfr},sr{pfr}(:,4)./(SDEN*abyv),'Color',[0 .45 .74])
plot(tr{pfr},sr{pfr}(:,5)./(SDEN*abyv),'k')
plot(tr{pfr},sr{pfr}(:,6)./(SDEN*abyv),'g')
plot(tr{pfr},sr{pfr}(:,10)./(SDEN*abyv),'m')
end
hold off
legend('N_{2*}','N_*','H_*','NH_{3*}','NH_{2*}','NH','\theta_*')
ss = sr{pfrnodes}(end,[1:6 10])/(SDEN*abyv);
fprintf('------------------------------------------------------------\n')
fprintf('                      Surface Species\n')
fprintf('------------------------------------------------------------\n')
fprintf('  N2       N        H2      NH3      NH2       NH       *\n')
fprintf('------   ------   ------   ------   ------   ------   ------\n')
fprintf('%6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f\n',ss)
fprintf('------------------------------------------------------------\n')