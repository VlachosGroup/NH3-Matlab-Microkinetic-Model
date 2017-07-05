function [ kf,kb,HORT,CpOR ] = amm_kinetics( T,T_gas,s )
% Ammonia Synthesis Microkinetic Model
%  Kinetic rate constants
%
% Input:
%   s      Species concentrations
%   T      Temperature, current (K)
%   T_orig Temperature, initial (k)
%
% Output:
%   kf  Forward reaction rate constants
%   kb  Backward reaction rate constants
%
% Species list:
%   s(1)    N2*
%   s(2)    N*
%   s(3)    H*
%   s(4)    NH3*
%   s(5)    NH2*
%   s(6)    NH*
%   s(7)    N2
%   s(8)    H2
%   s(9)    NH3
%   s(10)   * (Vacant site)
% 
%%
% Kinetic rate constants
global abyv T_ref beta Stoic_gas Stoic MWON A Stick MW_N2 MW_H2...
       MW_NH3 R_e R_k R Ea SDTOT

[Ea,A6_BEP,~] = amm_BEP_LSR(T,Stoic,Ea,s);
A6_Cov = amm_coverage(s);
kf=zeros(7,1);

kf(1) = 1*(Stick(1)/(1-MWON*Stick(1)/2))/SDTOT*((T_gas/T_ref)^beta(1))* ...
        sqrt(R_k*T_gas/(2*pi*MW_N2)) * exp(-Ea(1)/(R_e*T_gas));               % N2(gas)   +  Ru(Step) <--> N2(Step)
kf(2) = 1*A(1)*((T/T_ref)^beta(2))/abyv * exp(-Ea(2)/(R_e*T));                % N2(Step)  +  Ru(Step) <--> 2N(Step)
kf(3) = 1*(Stick(2)/(1-MWON*Stick(2)/2))/(abyv*SDTOT^2)*((T_gas/T_ref)^beta(3))* ...
        sqrt(R_k*T_gas/(2*pi*MW_H2)) * exp(-Ea(3)/(R_e*T_gas));               % H2(gas)   + 2Ru(Step) <--> 2H(Step)
kf(4) = 1*A(2)*((T/T_ref)^beta(4))/abyv * exp(-Ea(4)/(R_e*T));                % NH3(Step) +  Ru(Step) <--> NH2(Step) + H(Step)
kf(5) = 1*A(3)*((T/T_ref)^beta(5))/abyv * exp(-Ea(5)/(R_e*T));                % NH2(Step) +  Ru(Step) <--> NH(Step)  + H(Step)
kf(6) = 1*A(4)*((T/T_ref)^beta(6))/abyv * exp(-Ea(6)/(R_e*T));                % NH(Step)  +  Ru(Step) <--> N(Step)   + H(Step)
kf(7) = 1*(Stick(3)/(1-MWON*Stick(3)/2))/SDTOT*((T_gas/T_ref)^beta(7))* ...
        sqrt(R_k*T_gas/(2*pi*MW_NH3)) * exp(-Ea(7)/(R_e*T_gas));              % NH3(gas)  +  Ru(Step) <--> NH3(Step)
[CpOR,HORT,SOR,~] = amm_thermo(T,A6_BEP,A6_Cov);
E = SOR - HORT;
E_e = E * (Stoic)';
Kp = exp(E_e)';
Kc = Kp .* (1/(R*T)).^(sum(Stoic_gas,2));
kb = kf./Kc;
end

