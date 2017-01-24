function [ Ea,Zero_BE,omega ] = amm_BES_LSR( T,Stoic )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

R_e = 1.987e-3;            % Gas constant, (kcal/mol K)

%% Liner Scaling Relationships

%  Qi = Qi,ref + ALPHAi * (Q_target - Q_ref)

Q_ref = 102.35;  % N binding energy on reference metal [kcal/mol]

Q_target = 133.8; % N binding energy

% ALPHAi LSR slopes
alpha(1)  = 0.62036;  %  N2  [Terrace]
alpha(2)  = 1;        %  N   [Terrace]
alpha(3)  = 0.17;     %  H   [Terrace]
alpha(4)  = 0.14;     %  NH3 [Terrace]
alpha(5)  = 0.41;     %  NH2 [Terrace]
alpha(6)  = 0.71;     %  NH  [Terrace]
alpha(7)  = 0.62036;  %  N2  [Step]
alpha(8)  = 1.057;    %  N   [Step]
alpha(9)  = 0.18;     %  H   [Step]
alpha(10) = 0.14;     %  NH3 [Step]
alpha(11) = 0.391;    %  NH2 [Step]
alpha(12) = 0.708;    %  NH  [Step]

% Qi,ref (zero coverage reference binding energy of the species) [kcal/mol]
Qi_ref(1)  = -2.0779;  %  N2  [Terrace]
Qi_ref(2)  = Q_ref;    %  N   [Terrace]
Qi_ref(3)  = 57.4245;  %  H   [Terrace]
Qi_ref(4)  = 12.2999;  %  NH3 [Terrace]
Qi_ref(5)  = 45.8833;  %  NH2 [Terrace]
Qi_ref(6)  = 82.5372;  %  NH  [Terrace]
Qi_ref(7)  = 9.451;    %  N2  [Step]
Qi_ref(8)  = 106.224;  %  N   [Step]
Qi_ref(9)  = 58.0824;  %  H   [Step]
Qi_ref(10) = 22.6759;  %  NH3 [Step]
Qi_ref(11) = 63.9298;  %  NH2 [Step]
Qi_ref(12) = 91.8554;  %  NH  [Step]

Q = Qi_ref + alpha * (Q_target - Q_ref);
Zero_BE = (alpha * (Q_target - Q_ref))'/R_e;

%% Bronsted-Evans-Polanyi Relationships for activation barriers from Hrxn

%  (Ea)=m(deltaHrxn)+b

%  m coefficients
m(1) = 0.681;  %N2 dissociation (Terrace)
m(2) = 0.69;   %N2 dissociation (Step)
m(3) = 0.29;   %NH dehydrogenation
m(4) = 0.52;   %NH2 dehydrogenation
m(5) = 0.71;   %NH3 dehydrogenation

%  b constant

b(1) = 54.27;   %N2 dissociation (Terrace)
b(2) = 40.42;   %N2 dissociation (Step)
b(3) = 23.23;   %NH dehydrogenation
b(4) = 19.78;   %NH2 dehydrogenation
b(5) = 23.69;   %NH3 dehydrogenation

[HORT,SOR,GORT] = amm_thermo(T,Zero_BE);
HRXN = HORT * Stoic'*T*R_e;
Ea(1) = max(m(2) * HRXN(2) + b(2),-HRXN(2));
Ea(2) = max(m(5) * HRXN(4) + b(5),-HRXN(4));
Ea(3) = max(m(4) * HRXN(5) + b(4),-HRXN(5));
Ea(4) = max(m(3) * HRXN(6) + b(3),-HRXN(6));

omega_default = 0;
omega(1)  = omega_default;   % N2(gas)      +  Ru(Terrace) <--> N2(Terrace)
omega(2)  = m(1);            % N2(Terrace)  +  Ru(Terrace) <--> 2N(Terrace)
omega(3)  = omega_default;   % H2(Terrace)  + 2Ru(Terrace) <--> 2H(Terrace)
omega(4)  = m(3);            % NH3(Terrace) +  Ru(Terrace) <--> NH2(Terrace) + H(Terrace)
omega(5)  = m(4);            % NH2(Terrace) +  Ru(Terrace) <--> NH(Terrace)  + H(Terrace)
omega(6)  = m(5);            % NH(Terrace)  +  Ru(Terrace) <--> N(Terrace)   + H(Terrace)
omega(7)  = omega_default;   % NH3(gas)     +  Ru(Terrace) <--> NH3(Terrace)
omega(8)  = omega_default;   % N2(gas)   +  Ru(Step) <--> N2(Step)
omega(9)  = m(2);            % N2(Step)  +  Ru(Step) <--> 2N(Step)
omega(10) = omega_default;   % H2(Step)  + 2Ru(Step) <--> 2H(Step)
omega(11) = m(3);            % NH3(Step) +  Ru(Step) <--> NH2(Step) + H(Step)
omega(12) = m(4);            % NH2(Step) +  Ru(Step) <--> NH(Step)  + H(Step)
omega(13) = m(5);            % NH(Step)  +  Ru(Step) <--> N(Step)   + H(Step)
omega(14) = omega_default;   % NH3(gas)  +  Ru(Step) <--> NH3(Step)
omega = omega';

end

