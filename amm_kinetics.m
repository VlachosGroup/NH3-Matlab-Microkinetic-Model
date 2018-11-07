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
       MW_NH3 R_e R_k R Ea SDTOT strain

[Ea,A6_LSR,A6_Strain,~] = amm_BEP_LSR(T,Stoic,Ea,s,strain);
A6_Cov = amm_coverage(s);

kf=zeros(7,1);
kf(1) = 1*(Stick(1)/(1-MWON*Stick(1)/2))/SDTOT*((T_gas/T_ref)^beta(1))* ...
        sqrt(R_k*T_gas/(2*pi*MW_N2)) * exp(-Ea(1)/(R_e*T_gas));               % N2   +  * <--> N2*
kf(2) = 1*A(1)*((T/T_ref)^beta(2))/abyv * exp(-Ea(2)/(R_e*T));                % N2*  +  * <--> 2N*
kf(3) = 1*(Stick(2)/(1-MWON*Stick(2)/2))/(abyv*SDTOT^2)*((T_gas/T_ref)^beta(3))* ...
        sqrt(R_k*T_gas/(2*pi*MW_H2)) * exp(-Ea(3)/(R_e*T_gas));               % H2   + 2* <--> 2H*
kf(4) = 1*A(2)*((T/T_ref)^beta(4))/abyv * exp(-Ea(4)/(R_e*T));                % NH3* +  * <--> NH2* + H*
kf(5) = 1*A(3)*((T/T_ref)^beta(5))/abyv * exp(-Ea(5)/(R_e*T));                % NH2* +  * <--> NH*  + H*
kf(6) = 1*A(4)*((T/T_ref)^beta(6))/abyv * exp(-Ea(6)/(R_e*T));                % NH*  +  * <--> N*   + H*
kf(7) = 1*(Stick(3)/(1-MWON*Stick(3)/2))/SDTOT*((T_gas/T_ref)^beta(7))* ...
        sqrt(R_k*T_gas/(2*pi*MW_NH3)) * exp(-Ea(7)/(R_e*T_gas));              % NH3  +  * <--> NH3*
[CpOR,HORT,SOR,~] = amm_thermo(T,A6_LSR,A6_Cov,A6_Strain);
GORT = HORT - SOR;
GORT_e = GORT * (Stoic)';
Kp = exp(-GORT_e)';
Kc = Kp .* (1/(R*T)).^(sum(Stoic_gas,2));
kb = kf./Kc;
end

