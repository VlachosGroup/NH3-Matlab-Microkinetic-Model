function [ kf,kb ] = amm_kinetics( T,T_orig,s )
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
global V SDEN abyv T_ref beta omega Stoic_surf Stoic_gas Stoic MWON
R_e = 1.987e-3;            % Gas constant, (kcal/mol K)
R_k = 8.31451e7;             % Gas constant, (g cm2/mol K s)
kf = zeros(7,1);
kb = zeros(7,1);
Stick_N2 = 0.5;                % Sticking coefficient
Stick_H2 = 0.5;                % Sticking coefficient
Stick_NH3 = 0.5;               % Sticking coefficient
M_N2  = 28.0134;             % Molecular weight (g/mole)
M_H2  = 2.01594;             % Molecular weight (g/mole)
M_NH3 = 17.03061;            % Molecular weight (g/mole)
cat_surf = abyv*V;
%A = abyv*6.65e19*((T/T_ref)^beta);              % Pre-Exponential factor (1/s)
[Ea,Zero_BE,omega] = amm_BEP_LSR(T,Stoic);
k_cov = amm_coverage(s);
%k_cov =zeros(10,1);
%omega = ones(14,1)*0.5;
Ea_kf = (omega(8:14)).*(Stoic*k_cov);
kf=zeros(7,1);
kf(1) = abyv^2*(Stick_N2/(1-MWON*Stick_N2/2))*((T/T_ref)^beta)* ...
        sqrt(R_k*T/(2*pi*M_N2))* exp(0*Ea_kf(1)/(R_e*T));                   % N2(gas)   +  Ru(Step) <--> N2(Step)
kf(2) = abyv*1.16e18*((T/T_ref)^beta) * exp((-Ea(1) + Ea_kf(2))/(R_e*T));   % N2(Step)  +  Ru(Step) <--> 2N(Step)
kf(3) = abyv^3*(Stick_H2/(1-MWON*Stick_H2/2)/SDEN)*((T/T_ref)^beta)* ...
        sqrt(R_k*T/(2*pi*M_H2))* exp(0*Ea_kf(3)/(R_e*T));                   % H2(Step)  + 2Ru(Step) <--> 2H(Step)
kf(4) = abyv*2.05e19*((T/T_ref)^beta) * exp((-Ea(2) + Ea_kf(4))/(R_e*T));   % NH3(Step) +  Ru(Step) <--> NH2(Step) + H(Step)
kf(5) = abyv*1.06e20*((T/T_ref)^beta) * exp((-Ea(3) + Ea_kf(5))/(R_e*T));   % NH2(Step) +  Ru(Step) <--> NH(Step)  + H(Step)
kf(6) = abyv*8.38e19*((T/T_ref)^beta) * exp((-Ea(4) + Ea_kf(6))/(R_e*T));   % NH(Step)  +  Ru(Step) <--> N(Step)   + H(Step)
kf(7) = abyv^2*(Stick_NH3/(1-MWON*Stick_NH3/2))*((T/T_ref)^beta)* ...
        sqrt(R_k*T/(2*pi*M_NH3))* exp(0*Ea_kf(7)/(R_e*T));                  % NH3(gas)  +  Ru(Step) <--> NH3(Step)
[HORT,SOR,GORT] = amm_thermo(T,Zero_BE);
%E = SOR - (HORT - k_cov'/(R_e*T));
E = SOR - HORT;
%E_e = E * (Stoic)';
E_e = E * (Stoic)' + (Stoic*k_cov)'/(R_e*T);
Kp = exp(E_e)';
R = 82.057338; %(cm3 atm/K mole)
Kc = Kp .* (1/(R*T)).^(sum(Stoic_gas,2));
kb = kf./Kc;
end

