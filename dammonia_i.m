function [ d2s ] = dammonia( ~,s )
% Ammonia Synthesis Microkinetic Model
% Jacobian Matrix to assist stiff solvers
%
% Input:
%   s   Species concentrations
%   T   Temperature (K)
%
% Output:
%   d2s  Jacobian matrix
%
% Species list:
%   [s(1),ds(1)]    H2
%   [s(2),ds(2)]    N2
%   [s(3),ds(3)]    H*
%   [s(4),ds(4)]    N*
%   [s(5),ds(5)]    NH*
%   [s(6),ds(6)]    NH2*
%   [s(7),ds(7)]    NH3*
%   [s(8),ds(8)]    NH3
%   [s(9),ds(9)]    * (Vacant site)
% 
% Kinetic rate constants
%[kf,kb]=amm_kinetics(T,s);  % Obtain kinetics rate constants
%
% Reaction network
%
global kb;
global kf;
d2s      = zeros(10,10);
d2s = [-kb(1)-kf(2)*s(10),(2)*kb(2)*s(2),0,0,0,0,kf(1)*s(10),0,0,kf(1)*s(7)-kf(2)*s(1);...
    (2)*kf(2)*s(10),-(4)*kb(2)*s(2)-kb(6)*s(3),-kb(6)*s(2),0,0,kf(6)*s(10),0,0,0,(2)*kf(2)*s(1)+kf(6)*s(6);...
    0,-kb(6)*s(3),-(4)*kb(3)*s(3)-kb(6)*s(2)-kb(4)*s(5)-kb(5)*s(6),kf(4)*s(10),kf(5)*s(10)-kb(4)*s(3),kf(6)*s(10)-kb(5)*s(3),0,(2)*kf(3)*s(10)^(2),0,kf(4)*s(4)+kf(5)*s(5)+kf(6)*s(6)+(4)*kf(3)*s(8)*s(10);...
    0,0,kb(4)*s(5),-kb(7)-kf(4)*s(10),kb(4)*s(3),0,0,0,kf(7)*s(10),kf(7)*s(9)-kf(4)*s(4);...
    0,0,kb(5)*s(6)-kb(4)*s(5),kf(4)*s(10),-kb(4)*s(3)-kf(5)*s(10),kb(5)*s(3),0,0,0,kf(4)*s(4)-kf(5)*s(5);...
    0,kb(6)*s(3),kb(6)*s(2)-kb(5)*s(6),0,kf(5)*s(10),-kb(5)*s(3)-kf(6)*s(10),0,0,0,kf(5)*s(5)-kf(6)*s(6);...
    kb(1),0,0,0,0,0,-kf(1)*s(10),0,0,-kf(1)*s(7);...
    0,0,(2)*kb(3)*s(3),0,0,0,0,-kf(3)*s(10)^(2),0,-(2)*kf(3)*s(8)*s(10);...
    0,0,0,kb(7),0,0,0,0,-kf(7)*s(10),-kf(7)*s(9);...
    kb(1)-kf(2)*s(10),(2)*kb(2)*s(2)+kb(6)*s(3),(4)*kb(3)*s(3)+kb(6)*s(2)+kb(4)*s(5)+kb(5)*s(6),kb(7)-kf(4)*s(10),kb(4)*s(3)-kf(5)*s(10),kb(5)*s(3)-kf(6)*s(10),-kf(1)*s(10),-(2)*kf(3)*s(10)^(2),-kf(7)*s(10),-kf(2)*s(1)-kf(1)*s(7)-kf(4)*s(4)-kf(5)*s(5)-kf(6)*s(6)-kf(7)*s(9)-(4)*kf(3)*s(8)*s(10)];
end

