function [ ds ] = ammonia( t,s )
%%          -------------------------------------------------
%                        NH3  Micro-kinetic model
%                         Vlachos Research Group
%                 Chemical and Biomolecular Egineering
%                         University of Delaware
%
%             Gerhard R Wittreich, P.E.  (February 10, 2017)
%           --------------------------------------------------
%
%
%  ammonia.m      : Reaction ODE's and dynamical equations
%
%      requires: 
%                amm_kinetics.m : Provides forward and reverse rate constants
%
% Input:
%   s   Species concentrations
%   T   Temperature (K)
%
% Output:
%   ds  Reaction rate for each species
%
% Species list:
%   [s(1),ds(1)]    N2*
%   [s(2),ds(2)]    N*
%   [s(3),ds(3)]    H*
%   [s(4),ds(4)]    NH3*
%   [s(5),ds(5)]    NH2*
%   [s(6),ds(6)]    NH*
%   [s(7),ds(7)]    N2
%   [s(8),ds(8)]    H2
%   [s(9),ds(9)]    NH3
%   [s(10),ds(10)]   * (Vacant site)
%   [s(11),ds(11)]  Catalyst surface temperature
%   [s(12),ds(12)]  Gas temperature
%
global kf kb T_orig T_pulse T_gas V Q_in c_N2 c_H2 c_NH3 Isobaric...
       Moles_SiO2_Heated R_e Cp_SiO2_NIST pulse abyv q_constant q_pulse...
       T_func RR surf_cat strain_pulse strain tspan

T_gas = s(12);
T = s(11);
% h_cat = 2.4e-6; % Catalyst heat transfer coefficient [kcal/cm2 s K]
% switch pulse
%     case 0
%     q = q_constant;
%     T = s(11);
%     case 1
%     q = sin(pulstran(t/2-floor(t/2),[0:1:1],'tripuls',0.002/2).^2*pi/2)*q_pulse;
%     %T = T_func(t);
%     T = s(11);
% end
if strain_pulse
    period = 0.1; % Pulse period [sec]
    magnitude = 0.04; % Pulse amplitude
    strain = sin(2*pi*(t-tspan)/period)*magnitude
end
[kf,kb,HORT,CpOR]=amm_kinetics(T,T_gas,s);  % Obtain kinetics rate constants
[~,~,HORT_feed,~]=amm_kinetics(T_orig,T_orig,s);
[~,~,HORT_gas,~]=amm_kinetics(T_gas,T_orig,s);

% Adjust reactor outflow to maintain an isobaric reactor if Isobaric = 1
Q_r = (kb(1)*s(1)   - kf(1)*s(7)*s(10) +...
       kb(3)*s(3)^2 - kf(3)*s(8)*s(10)^2 +...
       kb(7)*s(4)   - kf(7)*s(9)*s(10))/sum(s(7:9))*V;
Q_out = Q_in + Q_r*Isobaric;
%
% Reaction network
%
RR(1,2) = kf(1)*s(7)*s(10);
RR(1,1) = kb(1)*s(1);
RR(1,3) = RR(1,1) - RR(1,2);
RR(2,2) = kf(2)*s(1)*s(10);
RR(2,1) = kb(2)*s(2)^2;
RR(2,3) = RR(2,1) - RR(2,2);
RR(3,2) = kf(3)*s(8)*s(10)^2;
RR(3,1) = kb(3)*s(3)^2;
RR(3,3) = RR(3,1) - RR(3,2);
RR(6,1) = kf(4)*s(4)*s(10);
RR(6,2) = kb(4)*s(5)*s(3);
RR(6,3) = RR(6,1) - RR(6,2);
RR(5,1) = kf(5)*s(5)*s(10);
RR(5,2) = kb(5)*s(6)*s(3);
RR(5,3) = RR(5,1) - RR(5,2);
RR(4,1) = kf(6)*s(6)*s(10);
RR(4,2) = kb(6)*s(2)*s(3);
RR(4,3) = RR(4,1) - RR(4,2);
RR(7,1) = kf(7)*s(9)*s(10);
RR(7,2) = kb(7)*s(4);
RR(7,3) = RR(7,1) - RR(7,2);
ds    = zeros(11,1);
ds(1) = kf(1)*s(7)*s(10)         - kb(1)*s(1) + ...
        kb(2)*s(2)^2             - kf(2)*s(1)*s(10);             % dN2*/dt
ds(2) = 2*kf(2)*s(1)*s(10)       - 2*kb(2)*s(2)^2 + ...
        kf(6)*s(6)*s(10)         - kb(6)*s(2)*s(3);              % dN*/dt
ds(3) = 2*kf(3)*s(8)*s(10)^2     - 2*kb(3)*s(3)^2 +...
        kf(4)*s(4)*s(10)         - kb(4)*s(5)*s(3) +...
        kf(5)*s(5)*s(10)         - kb(5)*s(6)*s(3) +...
        kf(6)*s(6)*s(10)         - kb(6)*s(2)*s(3);              % dH*/dt
ds(4) = kb(4)*s(5)*s(3)          - kf(4)*s(4)*s(10) +...
        kf(7)*s(9)*s(10)         - kb(7)*s(4);                   % dNH3*/dt
ds(5) = kf(4)*s(4)*s(10)         - kb(4)*s(5)*s(3) +...
        kb(5)*s(6)*s(3)          - kf(5)*s(5)*s(10);             % dNH2*/dt
ds(6) = kf(5)*s(5)*s(10)         - kb(5)*s(6)*s(3) +...
        kb(6)*s(2)*s(3)          - kf(6)*s(6)*s(10);             % dNH*/dt
ds(7) = kb(1)*s(1)               - kf(1)*s(7)*s(10) +...
        Q_in/V*c_N2 - Q_out/V*s(7);                              % dN2/dt
ds(8) = kb(3)*s(3)^2             - kf(3)*s(8)*s(10)^2 +...
        Q_in/V*c_H2 - Q_out/V*s(8);                              % dH2/dt
ds(9) = kb(7)*s(4)               - kf(7)*s(9)*s(10) +...
        Q_in/V*c_NH3 - Q_out/V*s(9);                             % dNH3/dt
ds(10)= kb(1)*s(1)               - kf(1)*s(7)*s(10) +...
        kb(2)*s(2)^2             - kf(2)*s(1)*s(10) +...
        2*kb(3)*s(3)^2           - 2*kf(3)*s(8)*s(10)^2 +...
        kb(4)*s(5)*s(3)          - kf(4)*s(4)*s(10) +...
        kb(5)*s(6)*s(3)          - kf(5)*s(5)*s(10) +...
        kb(6)*s(2)*s(3)          - kf(6)*s(6)*s(10) +...
        kb(7)*s(4)               - kf(7)*s(9)*s(10);             % d*/dt
ds(11)= 0;%(q + (HORT(1:6)*R_e*T)*ds(1:6)*V + ...
%          Q_in *R_e*T*(c_N2*HORT(7) + c_H2*HORT(8) + c_NH3*HORT(9)) -...
%          Q_out*R_e*T*(s(7)*HORT(7) + s(8)*HORT(8) + s(9)*HORT(9)) +...
%          h_cat*surf_cat*V*(T_gas - T))/...
%          (Moles_SiO2_Heated*(Cp_SiO2_NIST*...
%          [1 T/1000 (T/1000)^2 (T/1000)^3 1/(T/1000)^2]'));        % dT_cat/dt
ds(12)= 0;%(Q_in*R_e*T_orig*(c_N2*HORT_feed(7) + c_H2*HORT_feed(8) + c_NH3*HORT_feed(9))-...
%          Q_in*R_e*T*(c_N2*HORT(7) + c_H2*HORT(8) + c_NH3*HORT(9))+...
%          Q_out*R_e*T*(s(7)*HORT(7)    + s(8)*HORT(8) + s(9)*HORT(9)) -...
%          Q_out*R_e*T_gas*(s(7)*HORT_gas(7) + s(8)*HORT_gas(8)  + s(9)*HORT_gas(9)) -...
%          h_cat*surf_cat*V*(T_gas - T))/...
%          (V*R_e*(s(7)*CpOR(7)+s(8)*CpOR(8)+s(9)*CpOR(9)));       % dT_gas/dt
%           if t>6
%               ff=0;
%           end
end
