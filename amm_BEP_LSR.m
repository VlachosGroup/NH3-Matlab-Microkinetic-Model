function [ Ea,A6_LSR,A6_Strain,Q ] = amm_BEP_LSR( T,Stoic,Ea,s,strain)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

global R_e Q_name STYPE_TERRACE

%% Liner Scaling Relationships

%  Qi = Qi,ref + ALPHAi * (Q_target - Q_ref)

Q_ref = 102.35;  % N binding energy on reference metal (Pt) [kcal/mol]

%Q_target = 102.35; Q_name = 'Pt'; % N binding energy on Pt
%Q_target = 110.00; Q_name = 'Ni'; % N binding energy on Ni
%Q_target = 112.07; Q_name = 'Rh'; % N binding energy on Rh
%Q_target = 115.30; Q_name = 'Co'; % N binding energy on Co
Q_target = 134.21; Q_name = 'Ru'; % N binding energy on Ru
%Q_target = 136.75; Q_name = 'Fe'; % N binding energy on Fe
%Q_target = 138.36; Q_name = 'Re'; % N binding energy on Re
%Q_target = 154.18; Q_name = 'Mo'; % N binding energy on Mo
%Q_target = 134.6595; Q_name = 'Unk';

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
A6_LSR = (alpha * (Q_target - Q_ref))'/R_e;

%% Catalyst surface strain

StrainCoef = [-0.04 -0.02 0.0 0.02 0.04];
StrainN2  = [0.15 0.11 0.00 -0.13 -0.21]*23.05875998694585/R_e;
StrainN   = [0.53 0.27 0.00 -0.30 -0.48]*23.05875998694585/R_e;
StrainH   = [0.12 0.05 0.00 -0.07 -0.10]*23.05875998694585/R_e;
StrainNH3 = [0.08 0.05 0.00 -0.11 -0.18]*23.05875998694585/R_e;
StrainNH2 = [0.12 0.07 0.00 -0.09 -0.15]*23.05875998694585/R_e;
StrainNH  = [0.26 0.14 0.00 -0.17 -0.27]*23.05875998694585/R_e;
A6_Strain = zeros(12,1);
A6_Strain(1)  = polyval(polyfit(StrainCoef,StrainN2,1), strain);
A6_Strain(2)  = polyval(polyfit(StrainCoef,StrainN,1), strain);
A6_Strain(3)  = polyval(polyfit(StrainCoef,StrainH,1), strain);
A6_Strain(4)  = polyval(polyfit(StrainCoef,StrainNH3,1), strain);
A6_Strain(5)  = polyval(polyfit(StrainCoef,StrainNH2,1), strain);
A6_Strain(6)  = polyval(polyfit(StrainCoef,StrainNH,1), strain);
A6_Strain(7:12) = A6_Strain(1:6);
%
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

A6_Cov = amm_coverage(s);
[~,HORT,~,~] = amm_thermo(T,A6_LSR,A6_Cov,A6_Strain);
HRXN = HORT * Stoic'*T*R_e;
Ea(2) = m(2) * HRXN(2) + b(2);
Ea(4) = m(5) * HRXN(4) + b(5);
Ea(5) = m(4) * HRXN(5) + b(4);
Ea(6) = m(3) * HRXN(6) + b(3);
end

