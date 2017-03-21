%%          -------------------------------------------------
%                        NH3  Micro-kinetic model
%                         Vlachos Research Group
%                 Chemical and Biomolecular Egineering
%                         University of Delaware
%
%             Gerhard R Wittreich, P.E.  (February 10, 2017)
%           --------------------------------------------------
%
%  Main program: 
%      requires: ammonia.m      : Reaction ODE's and dynamical equations
%                amm_kinetics.m : Provides forward and reverse rate constants
%                amm_thermo.m   : NASA polynomials provides enthalpy and entropy
%                amm_coverage.m : Calculates lateral interactions of surface species
%                amm_BEP_LSR.m  : Adjusts thermo data for target metal
%                                 catalyst, provides Ea and omega for each reaction
%%
clear
datetime
%  Set key model parameters
%
global T T_ref T_pulse T_orig beta P V Q_in c_N2 c_H2 c_NH3 abyv ...
       c_tot Stoic_surf Stoic_gas Stoic MWON Isobaric Ea A Stick R_e ...
       R_k R MW_N2 MW_H2 MW_NH3 SDEN_1 SDEN_2 SDTOT Moles_SiO2_Heated...
       Cp_SiO2_NIST pulse
T_orig = 700;                   % Reactor bulk temperature [K]
T = T_orig;                     % Initial reactor temperature [K]
P = 1;                          % Reactor pressure [atm]
beta = [0 1 0 1 1 1 0]';
Ea   = [0 0 0 0 0 0 0]';
A    = [1.16e18 2.05e19 1.06e20 8.38e19]';
Stick= [0.5 0.5 0.5]';
MWON = 0;                       % Motz-Wise Correction: 0 = Off (S); 1 = On (S/(1-S/2))
Isobaric = 1;                   % Isobaric Reactor (0=Isochoric, 1=Isobaric)
T_ref = 1;                      % Reference Temperature [K]
MW_H = 1.00797;                 %---
MW_N = 14.0067;                 %
MW_H2 = 2*MW_H;                 % Molecular mass [g/mol]
MW_N2 = 2*MW_N;                 %
MW_NH3 = MW_N + 3*MW_H;         %---
X_H2  = 0;                      %
X_N2  = 0;                      % Initial mole fraction in reactor feed
X_NH3 = 1;                      %
Y_H2  = X_H2 /(X_H2+X_N2+X_NH3);% Mole fractions
Y_N2  = X_N2 /(X_H2+X_N2+X_NH3);% normalized
Y_NH3 = X_NH3/(X_H2+X_N2+X_NH3);% to 1
abyv = 1500;                    % Catalyst loading (cm2 catalyst/cm3 reac volume)
V = 10;                         % Reactor volume (cm3)
Q_in = 10;                      % 0 = Batch Reactor,  Any other value = CSTR [cm3/s]
eps = 0.64;                     % Sphere packed volume
Cat_Rad = 3/abyv;               %
n_Cat = V*3/(4*pi*Cat_Rad^3);   %
V_Cat_Heat = 4/3*pi*(Cat_Rad^3-(0.9*Cat_Rad)^3)*n_Cat;
rho_Cat = 2.68;                 % Density SiO2 catalyst substrate (gm/cm3)
MW_SiO2 = 60.08;                % Molecular weight SiO2 [gm/mol]
Moles_SiO2_Heated = V_Cat_Heat*rho_Cat/MW_SiO2;
Cp_SiO2_NIST = [-1.452341 60.15189 -77.6282 40.2869 0.000609]/1000;
SDEN_1 = 2.1671e-09;            % Catalyst terrace site density (moles/cm2)
SDEN_2 = 4.4385e-10;            % Catalyst step site density (moles/cm2)
SDTOT = SDEN_1 + SDEN_2;        % Total catalyst site density (moles/cm2)
R_e = 1.987e-3;                 % Gas constant, (kcal/mol K)
R_k = 8.31451e7;                % Gas constant, (g cm2/mol K s)
R = 82.057;                     % Gas constant, (cm3 atm/K mol)
c_tot = P/(R*T);                % Total starting moles in reactor [mol/cm3]
c_H2 = Y_H2*c_tot;              % Moles H2 in reactor [mol/cm3]
c_N2 = Y_N2*c_tot;              % Moles N2 in reactor [mol/cm3]
c_NH3 = Y_NH3*c_tot;            % Moles NH3 in reactor [mol/cm3]
M_H2 = c_H2*V*MW_H2;            % Mass H2 in reactor [kg]
M_N2 = c_N2*V*MW_N2;            % Mass N2 in reactor [kg]
M_NH3 = c_NH3*V*MW_NH3;         % Mass NH3 in reactor [kg]
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
options0 = odeset ('MaxStep',0.003,'NonNegative',[1 2 3 4 5 6 7 8 9 10],...
                   'BDF','on','InitialStep',1e-10,'Stats','off',...
                   'AbsTol',1e-14,'RelTol',1e-12);
options1 = odeset ('MaxStep',0.0003,'NonNegative',[1 2 3 4 5 6 7 8 9 10],...
                   'BDF','on','InitialStep',1e-10,'Stats','off',...
                   'AbsTol',1e-14,'RelTol',1e-12);
options2 = odeset ('NonNegative',[1 2 3 4 5 6 7 8 9 10],'InitialStep',1e-10,...
                   'BDF','on','Stats','off','AbsTol',1e-14,'RelTol',1e-12);
tic;
s0 = [0 0 0 0 0 0 c_N2 c_H2 c_NH3 SDEN_2*abyv T T]; % Initial species concentrations
if ne(0,1)
    %T_pulse = T_orig;
    pulse = 0;
    tspan = 1200;%max(floor(5*V/Q_in),5);
    [t,s] = ode15s(@ammonia,[0 tspan],s0,options2);
    s(:,10) = (SDEN_2*abyv) - sum(s(:,1:6),2);
    s0 = s(end,1:12);
    clear t s
end
tstart = 0;
pfrnodes = 1;           % PFR capability is not implemented.  Must be 1.
for pfr=1:pfrnodes
    %T=700;
    %T_pulse = 800;
    pulse = 1;
    tspan = 1200;%max(floor(3*V/Q_in)+.5,2);
    [t,s] = ode15s(@ammonia,[tstart tspan+tstart],s0,options0);
    s(:,10) = (SDEN_2*abyv) - sum(s(:,1:6),2);
    tr{pfr}=t;
    sr{pfr}=s;
end
toc;
save('ammonia_temp.mat')
figure(1)
hold on
for pfr=1:pfrnodes
plot(tr{pfr},sr{pfr}(:,7)./sum(sr{pfr}(:,7:9),2),'b')
plot(tr{pfr},sr{pfr}(:,8)./sum(sr{pfr}(:,7:9),2),'r')
plot(tr{pfr},sr{pfr}(:,9)./sum(sr{pfr}(:,7:9),2),'g')
end
hold off
xlim([0 tspan])
ylim([0 1])
xlabel('Time [sec]')
ylabel('Mole fraction [gas]')
legend('N_2','H_2','NH_3')
NH3_Conv = 1-sr{1}(end,9)/c_NH3;
fprintf('\n------------------------------\n')
fprintf('NH3 Conversion = %6.2f\n',NH3_Conv*100)
fs = sr{pfrnodes}(end,7:9)./sum(sr{pfrnodes}(end,7:9),2);
fprintf('------------------------------\n')
fprintf('         Gas Species\n')
fprintf('------------------------------\n')
fprintf('   N2         H2         NH3\n')
fprintf('--------   --------   --------\n')
fprintf('%8.6f   %8.6f   %8.6f\n',fs)
figure(2)
hold on
for pfr=1:pfrnodes
plot(tr{pfr},sr{pfr}(:,1) ./(SDEN_2*abyv),'-b')
plot(tr{pfr},sr{pfr}(:,2) ./(SDEN_2*abyv),'-r')
plot(tr{pfr},sr{pfr}(:,3) ./(SDEN_2*abyv),'-c')
plot(tr{pfr},sr{pfr}(:,4) ./(SDEN_2*abyv),'-','Color',[0 .45 .74])
plot(tr{pfr},sr{pfr}(:,5) ./(SDEN_2*abyv),'-k')
plot(tr{pfr},sr{pfr}(:,6) ./(SDEN_2*abyv),'-g')
plot(tr{pfr},sr{pfr}(:,10)./(SDEN_2*abyv),'-m')
end
hold off
xlim([0 tspan])
ylim([0 1])
xlabel('Time [sec]')
ylabel('Surface coverage')
legend('N_{2*}','N_*','H_*','NH_{3*}','NH_{2*}','NH','\theta_*')
ss = sr{pfrnodes}(end,[1:6 10])/(SDEN_2*abyv);
fprintf('--------------------------------------------------------------------------\n')
fprintf('                             Surface Species\n')
fprintf('--------------------------------------------------------------------------\n')
fprintf('   N2         N          H2        NH3        NH2         NH         *\n')
fprintf('--------   --------   --------   --------   --------   --------   --------\n')
fprintf('%8.6f   %8.6f   %8.6f   %8.6f   %8.6f   %8.6f   %8.6f\n',ss)
fprintf('--------------------------------------------------------------------------\n')
hold off
figure(3)
hold on
for pfr=1:pfrnodes
plot(tr{pfr},sr{pfr}(:,11),'-b')
end
title('Catalyst Surface Temperature')
xlim([0 tspan])
xlabel('Time [sec]')
ylabel('Temperature [K]')
hold off
figure(4)
hold on
for pfr=1:pfrnodes
plot(tr{pfr},sr{pfr}(:,12),'-b')
end
title('Gas Temperature')
xlim([0 tspan])
xlabel('Time [sec]')
ylabel('Temperature [K]')
hold off