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
global SDEN abyv T_ref beta omega Stoic_gas Stoic MWON A Stick MW_N2 MW_H2...
       MW_NH3 R_e R_k R Ea
kf = zeros(7,1);
kb = zeros(7,1);
[Ea,Zero_BE,omega] = amm_BEP_LSR(T,Stoic,Ea);
k_cov = amm_coverage(s);
%k_cov =zeros(10,1);
Ea_kf = (omega(8:14)).*(Stoic*k_cov);
ff=1;
kf=zeros(7,1);
kf(1) = (Stick(1)/(1-MWON*Stick(1)/2))/SDEN*((T/T_ref)^beta(1))* ...
        sqrt(R_k*T/(2*pi*MW_N2))* exp((-Ea(1) + Ea_kf(1))/(R_e*T));           % N2(gas)   +  Ru(Step) <--> N2(Step)
kf(2) = A(1)*((T/T_ref)^beta(2))/abyv * exp((-Ea(2) + Ea_kf(2))/(R_e*T));     % N2(Step)  +  Ru(Step) <--> 2N(Step)
kf(3) = (Stick(2)/(1-MWON*Stick(2)/2))/(abyv*SDEN^2)*((T/T_ref)^beta(3))* ...
        sqrt(R_k*T/(2*pi*MW_H2))* exp((-Ea(3) + Ea_kf(3))/(R_e*T));           % H2(gas)   + 2Ru(Step) <--> 2H(Step)
kf(4) = A(2)*((T/T_ref)^beta(4))/abyv * exp((-Ea(4) + Ea_kf(4))/(R_e*T));     % NH3(Step) +  Ru(Step) <--> NH2(Step) + H(Step)
kf(5) = A(3)*((T/T_ref)^beta(5))/abyv * exp((-Ea(5) + Ea_kf(5))/(R_e*T));     % NH2(Step) +  Ru(Step) <--> NH(Step)  + H(Step)
kf(6) = A(4)*((T/T_ref)^beta(6))/abyv * exp((-Ea(6) + Ea_kf(6))/(R_e*T));     % NH(Step)  +  Ru(Step) <--> N(Step)   + H(Step)
kf(7) = (Stick(3)/(1-MWON*Stick(3)/2))/SDEN*((T/T_ref)^beta(7))* ...
        sqrt(R_k*T/(2*pi*MW_NH3))* exp((-Ea(7) + Ea_kf(7))/(R_e*T));          % NH3(gas)  +  Ru(Step) <--> NH3(Step)
[HORT,SOR,GORT] = amm_thermo(T,Zero_BE);
E = SOR - HORT;
E_e = E * (Stoic)' + (Stoic*k_cov)'/(R_e*T);
Kp = exp(E_e)';
Kc = Kp .* (1/(R*T)).^(sum(Stoic_gas,2));
kb = kf./Kc;
end

