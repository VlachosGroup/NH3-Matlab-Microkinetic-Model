function [ k_cov ] = amm_coverage( s )
% Ammonia Synthesis Microkinetic Model
%  Kinetic rate constants
%
% Input:
%   s   Species concentrations
%
% Output:
%   k_cov  Catalyst coverage effects
%
% Species list:
%     s(1)             N2*
%     s(2)             N*
%     s(3),k_cov(1)    H*
%     s(4),k_cov(2)    NH3*
%     s(5),k_cov(3)    NH2*
%     s(6),k_cov(4)    NH*
%     s(7),k_cov(5)    N2
%     s(8)             H2
%     s(9)             NH3
%     s(10)             * (Vacant site)
%
global SDEN
%            N*        H*         NH3*       NH2*        NH*
Effects = zeros(10);
Effects(2:6,2:6) = [-46.9	-17.71  	-25.1       -20.71      -48.66;...
           -17.71	-6.6875 	-9.4781 	-7.8203 	-18.3746;...
           -25.1	-9.4781     -13.433     -11.0836	-26.0419;...
           -20.71	-7.8203     -11.0836	-9.1451 	-21.4872;...
           -48.66	-18.3746	-26.0419   	-21.4872	-50.486];
           
k_cov = (Effects'*(s/SDEN));
end

