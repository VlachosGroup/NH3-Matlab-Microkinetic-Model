function [tr, sr, RR] = amm_main(T_in, Graph)
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
%clear
datetime,T_in
%  Set key model parameters
%
global T T_ref T_pulse T_orig beta P V Q_in c_N2 c_H2 c_NH3 abyv ...
    c_tot Stoic_surf Stoic_gas Stoic MWON Isobaric Ea A Stick R_e ...
    R_k R MW_N2 MW_H2 MW_NH3 SDEN_1 SDEN_2 SDTOT Moles_SiO2_Heated...
    Cp_SiO2_NIST pulse q_constant q_pulse T_func RR Q_name
T_orig = T_in;                   % Reactor bulk temperature [K]
T = T_orig;                     % Initial catalyst temperature [K]
T_gas = T_orig;                 % Initial gas temperature [K]
P = 1.0;                        % Reactor pressure [atm]
q_constant = 0.0015;
q_pulse = 0.281636169371400;
beta = [0 1 0 1 1 1 0]';
Ea   = [0 0 0 0 0 0 0]';
A    = [1.16e18 2.05e19 1.06e20 8.38e19]';
Stick= [0.5 0.5 0.5]';
MWON = 0;                       % Motz-Wise Correction: 0 = Off (S); 1 = On (S/(1-S/2))
Isobaric = 1;                   % Isobaric Reactor (0=Isochoric, 1=Isobaric)
T_ref = 1;                      % Reference Temperature [K]
MW_H = 1.00797;                 % ---
MW_N = 14.0067;                 %
MW_H2 = 2*MW_H;                 % Molecular mass [g/mol]
MW_N2 = 2*MW_N;                 %
MW_NH3 = MW_N + 3*MW_H;         % ---
X_H2  = 0;                      %
X_N2  = 0;                      % Initial mole fraction in reactor feed
X_NH3 = 1;                      %
Y_H2  = X_H2 /(X_H2+X_N2+X_NH3);% Mole fractions
Y_N2  = X_N2 /(X_H2+X_N2+X_NH3);% normalized
Y_NH3 = X_NH3/(X_H2+X_N2+X_NH3);% to 1
abyv = 1500;                    % Catalyst loading (cm2 catalyst/cm3 reac volume)
V = 1;                        % Reactor volume (cm3)
Q_in = 145/60*T_orig/298.15/P;  % 0 = Batch Reactor,  Any other value = CSTR [cm3/s]
Q_in = 10;                       % 0 = Batch Reactor,  Any other value = CSTR [cm3/s]
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
options0 = odeset ('MaxStep',0.003,'NonNegative',[1 2 3 4 5 6 7 8 9 10 11 12],...
    'BDF','on','InitialStep',1e-10,'Stats','off',...
    'AbsTol',1e-14,'RelTol',1e-12);
options1 = odeset ('MaxStep',0.0001,'NonNegative',[1 2 3 4 5 6 7 8 9 10 11 12],...
    'BDF','on','InitialStep',1e-10,'Stats','off',...
    'AbsTol',1e-14,'RelTol',1e-12);
options2 = odeset ('NonNegative',[1 2 3 4 5 6 7 8 9 10 11 12],'InitialStep',1e-10,...
    'BDF','on','Stats','off','AbsTol',1e-14,'RelTol',1e-12);
tic;
s0 = [0 0 0 0 0 0 c_N2 c_H2 c_NH3 SDEN_2*abyv T T_gas]; % Initial species concentrations
if ne(0,1)
    T_pulse = T_orig;
    pulse = 0;
    tspan = max(floor(5*V/Q_in),5);
    sol = ode15s(@ammonia,[0 tspan],s0,options2);
    s0 = sol.y(:,end)';
end
tstart = 0;
pfrnodes = 1;           % PFR capability is not implemented.  Must be 1.
for pfr=1:pfrnodes
    %T=700;
    T_pulse = T_orig;
    pulse = 0;
    tspan2 = max(floor(3*V/Q_in)+.5,2);
    t = [tspan:0.0001:tspan2+tspan+1];
    T_func = griddedInterpolant(t,sin(pulstran(t/2-floor(t/2),[0:0.1:1],'tripuls',0.008).^3*pi/2)*(T_pulse-T_orig)+T_orig);
    sol2 = odextend(sol,@ammonia,tspan+tspan2,s0,options2);
    %s(:,10) = (SDEN_2*abyv) - sum(s(:,1:6),2);
    tr{pfr}=sol2.x';
    sr{pfr}=sol2.y';
end
toc;
if Graph
save('ammonia_temp_Ru_2sec.mat')
figure(1)
hold on
for pfr=1:pfrnodes
    plot(tr{pfr},sr{pfr}(:,7)./sum(sr{pfr}(:,7:9),2),'b')
    plot(tr{pfr},sr{pfr}(:,8)./sum(sr{pfr}(:,7:9),2),'r')
    plot(tr{pfr},sr{pfr}(:,9)./sum(sr{pfr}(:,7:9),2),'g')
end
hold off
xlim([0 tspan+tspan2])
ylim([0 1])
xlabel('Time [sec]')
ylabel('Mole fraction [gas]')
legend('N_2','H_2','NH_3')

switch pulse
    case 0
        Energy = q_constant;
        NH3_Conv = 1 - sr{1}(end,9)/c_NH3;
        Tf_cat = sr{1}(end,11);
        Tf_gas = sr{1}(end,12);
    case 1
        Energy = integral(@(t) sin(pulstran(t-floor(t),[0:1:1],'tripuls',0.01).^2*pi/2)*q_pulse,0,1000);
        NH3_Conv = 1-trapz(tr{1}(find(tr{1}==tspan):end),sr{1}(find(tr{1}==tspan):end,9))/...
            (tr{1}(end)-tr{1}(find(tr{1}==tspan)))/c_NH3;
        Tf_cat = trapz(tr{1}(find(tr{1}==tspan):end),sr{1}(find(tr{1}==tspan):end,11))/...
            (tr{1}(end)-tr{1}(find(tr{1}==tspan)));
        Tf_gas = trapz(tr{1}(find(tr{1}==tspan):end),sr{1}(find(tr{1}==tspan):end,12))/...
            (tr{1}(end)-tr{1}(find(tr{1}==tspan)));
end
H_conv = Energy/(Q_in*c_NH3*NH3_Conv);
fprintf('\n----------------------------------------\n')
fprintf('NH3 Conversion = %6.4f [%%]\n',NH3_Conv*100)
fprintf('Conv Enthalpy  = %6.4f [kcal/mol]\n',H_conv)
fprintf('Temperature (Feed)      = %7.4f [K]\n',T_orig)
fprintf('Temperature (Catalyst)  = %7.4f [K]\n',Tf_cat)
fprintf('Temperature (Gas)       = %7.4f [K]\n',Tf_gas)
fs = sr{pfrnodes}(end,7:9)./sum(sr{pfrnodes}(end,7:9),2);
fprintf('----------------------------------------\n')
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
xlim([0 tspan+tspan2])
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
    plot(tr{pfr},sr{pfr}(:,11),'r')
end
title('Catalyst Surface Temperature')
xlim([0 tspan+tspan2])
xlabel('Time [sec]')
ylabel('Temperature [K]')
hold off
figure(4)
hold on
for pfr=1:pfrnodes
    plot(tr{pfr},sr{pfr}(:,12),'r')
end
title('Gas Temperature')
xlim([0 tspan+tspan2])
xlabel('Time [sec]')
ylabel('Temperature [K]')
hold off
figure(5)
TOF = RR ./abyv ./(SDEN_1 + SDEN_2);
BG = barh(TOF);
set(BG(3),'DisplayName','Net','FaceColor','r');
set(BG(2),'DisplayName','Reverse','FaceColor','y');
set(BG(1),'DisplayName','Forward','FaceColor','b');
set(gca,'Xscale','log', 'XMinorTick', 'off')
set(gca,'yticklabel',{'N_2* \leftrightarrow N_2 + *';...
                      '2N* \leftrightarrow N_2* + *';...
                      '2H* \leftrightarrow H_2 + 2*';...
                      'NH* + * \leftrightarrow N* + H*';...
                      'NH_2* + * \leftrightarrow NH* + H*';...
                      'NH_3* + * \leftrightarrow NH_2* + H*';...
                      'NH_3 + * \leftrightarrow NH_3*'})
xlabel('Turnover Frequency (TOF) [s^{-1}]')
ylabel('Reaction Step')
title({'Ammonia Decomposition', ['Forward and Reverse Reaction Rates at ' ...
        num2str(sr{1}(end,12)) ' [K] on ' Q_name],['V_{Reactor} = ' ...
        num2str(V) ' cm^3     Q_{Feed} = ' num2str(Q_in) ...
        ' cm^3/s     \tau_{Reactor} = ' num2str(V/Q_in) ' seconds']})
legend('Forward', 'Reverse', 'Net', 'Location', 'best')
figure(6)
PEI = RR(:,1)./(RR(:,1)+RR(:,2));
plot(PEI,[1:length(PEI)],'o', 'MarkerFacecolor','b')
hold on
h=fill([0.45 0.45 0.55 0.55],[1 7 7 1],'y');
set(h,'facealpha',.1,'linestyle','none');
xlim([0,1])
set(gca,'yticklabel',{'N_2* \leftrightarrow N_2 + *';...
'2N* \leftrightarrow N_2* + *';...
'2H* \leftrightarrow H_2 + 2*';...
'NH* + * \leftrightarrow N* + H*';...
'NH_2* + * \leftrightarrow NH* + H*';...
'NH_3* + * \leftrightarrow NH_2* + H*';...
'NH_3 + * \leftrightarrow NH_3*'})
plot([0.45 0.45],[1,7],'--b')
plot([0.55 0.55],[1,7],'--b')
ylabel('Reaction Step')
xlabel('Partial Equilibrium Index')
hold off
end
end