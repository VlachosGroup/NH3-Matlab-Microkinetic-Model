function [ kf,kb ] = amm_kinetics( T,T_orig,P,s,Stoic,Stoic_surf,Stoic_gas )
% Ammonia Synthesis Microkinetic Model
%  Kinetic rate constants
%
% Input:
%   s   Species concentrations
%   T   Temperature (K)
%
% Output:
%   kf  Forward reaction rate constants
%   kb  Backward reaction rate constants
%
% Species list:
%     s(1)    H2
%     s(2)    N2
%     s(3)    H*
%     s(4)    N*
%     s(5)    NH*
%     s(6)    NH2*
%     s(7)    NH3*
%     s(8)    NH3
%     s(9)    * (Vacant site)
%     s(10)   N2*
% 

% Kinetic rate constants
global V SDEN abyv T_ref beta
R_e = 1.987e-3;            % Gas constant, (kcal/mol K)
R_k = 8.31451e7;             % Gas constant, (g cm2/mol K s)
kf = zeros(7,1);
kb = zeros(7,1);
omega = 0;
Stick_N2 = 1;                % Sticking coefficient
Stick_H2 = 1;                % Sticking coefficient
Stick_NH3 = 1;               % Sticking coefficient
M_N2  = 28.0134;             % Molecular weight (g/mole)
M_H2  = 2.01594;             % Molecular weight (g/mole)
M_NH3 = 17.03061;            % Molecular weight (g/mole)
cat_surf = abyv*V;
A = abyv*6.65e19*((T/T_ref)^beta);              % Pre-Exponential factor (1/s)

k_cov = amm_coverage(s);
%k_cov =zeros(10,1);
Ea_kf = omega*Stoic_surf*k_cov;
kf=zeros(7,1);
kf(1) = abyv*(Stick_N2/(1-Stick_N2/2))*((T_orig/T_ref)^beta)* sqrt(R_k*T_orig/(2*pi*M_N2))* exp(Ea_kf(1)/(R_e*T_orig));
kf(2) = A * exp((-36.7 + Ea_kf(2))/(R_e*T));
kf(3) = abyv^2*(Stick_H2/(1-Stick_H2/2)/SDEN)*((T_orig/T_ref)^beta)* sqrt(R_k*T_orig/(2*pi*M_H2))* exp(Ea_kf(3)/(R_e*T_orig));
kf(4) = A * exp((-21.9 + Ea_kf(4))/(R_e*T));
kf(5) = A * exp((-18.2 + Ea_kf(5))/(R_e*T));
kf(6) = A * exp((-21.4 + Ea_kf(6))/(R_e*T));
kf(7) = abyv*(Stick_NH3/(1-Stick_NH3/2))*((T_orig/T_ref)^beta)* sqrt(R_k*T_orig/(2*pi*M_NH3))* exp(Ea_kf(7)/(R_e*T_orig));
Kc = amm_thermo(T,P,Stoic,Stoic_gas,k_cov,omega);
%Kc(1)=1000;Kc(2)=1000;Kc(3)=1000;Kc(4)=0.00001;Kc(5)=0.00001;Kc(6)=0.00001;Kc(7)=0.001;
kb = kf./Kc;
end

