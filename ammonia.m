function [ ds ] = ammonia( t,y )
% Ammonia Synthesis Microkinetic Model
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
%
global kf kb T T_orig T_pulse V Q_in c_N2 c_H2 c_NH3 SDEN Isobaric abyv R

%y(10) = SDEN*abyv - sum(y(1:6));   % Determine vacancy site concentration

%T = sin(pulstran(t-floor(t/10)*10,[0:1:10],'tripuls',0.01).^2*pi/2)*(T_pulse-T_orig)+T_orig;


[kf,kb]=amm_kinetics(T,T_orig,y);  % Obtain kinetics rate constants
% Adjust reactor outflow to maintain an isobaric reactor if Isobaric = 1
Q_r = (kb(1)*y(1)   - kf(1)*y(7)*y(10) +...
       kb(3)*y(3)^2 - kf(3)*y(8)*y(10)^2 +...
       kb(7)*y(4)   - kf(7)*y(9)*y(10))/sum(y(7:9))*V;
Q_out = Q_in + Q_r*Isobaric;
%
% Reaction network
%
ds    = zeros(9,1);
ds(1) = kf(1)*y(7)*y(10)         - kb(1)*y(1) + ...
        kb(2)*y(2)^2             - kf(2)*y(1)*y(10);             % dN2*/dt
ds(2) = 2*kf(2)*y(1)*y(10)       - 2*kb(2)*y(2)^2 + ...
        kf(6)*y(6)*y(10)         - kb(6)*y(2)*y(3);              % dN*/dt
ds(3) = 2*kf(3)*y(8)*y(10)^2     - 2*kb(3)*y(3)^2 +...
        kf(4)*y(4)*y(10)         - kb(4)*y(5)*y(3) +...
        kf(5)*y(5)*y(10)         - kb(5)*y(6)*y(3) +...
        kf(6)*y(6)*y(10)         - kb(6)*y(2)*y(3);              % dH*/dt
ds(4) = kb(4)*y(5)*y(3)          - kf(4)*y(4)*y(10) +...
        kf(7)*y(9)*y(10)         - kb(7)*y(4);                   % dNH3*/dt
ds(5) = kf(4)*y(4)*y(10)         - kb(4)*y(5)*y(3) +...
        kb(5)*y(6)*y(3)          - kf(5)*y(5)*y(10);             % dNH2*/dt
ds(6) = kf(5)*y(5)*y(10)         - kb(5)*y(6)*y(3) +...
        kb(6)*y(2)*y(3)          - kf(6)*y(6)*y(10);             % dNH*/dt
ds(7) = kb(1)*y(1)               - kf(1)*y(7)*y(10) +...
        Q_in/V*c_N2 - Q_out/V*y(7);                              % dN2/dt
ds(8) = kb(3)*y(3)^2             - kf(3)*y(8)*y(10)^2 +...
        Q_in/V*c_H2 - Q_out/V*y(8);                              % dH2/dt
ds(9) = kb(7)*y(4)               - kf(7)*y(9)*y(10) +...
        Q_in/V*c_NH3 - Q_out/V*y(9);                             % dNH3/dt
ds(10)= kb(1)*y(1)               - kf(1)*y(7)*y(10) +...
        kb(2)*y(2)^2             - kf(2)*y(1)*y(10) +...
        2*kb(3)*y(3)^2           - 2*kf(3)*y(8)*y(10)^2 +...
        kb(4)*y(5)*y(3)          - kf(4)*y(4)*y(10) +...
        kb(5)*y(6)*y(3)          - kf(5)*y(5)*y(10) +...
        kb(6)*y(2)*y(3)          - kf(6)*y(6)*y(10) +...
        kb(7)*y(4)               - kf(7)*y(9)*y(10);             % d*/dt
end
