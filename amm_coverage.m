function [ A6_Cov ] = amm_coverage( s )
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
global SDEN_2 abyv R_e
Effects = zeros(10);
%                      N*          H*         NH3*        NH2*         NH*
Effects(2:6,2:6) = [-47.0179	-17.7545  	-25.1631    -20.7620    -48.7823 ;...
                    -17.7545	-6.7043 	-9.5019 	-7.8400     -18.4208 ;...
                    -25.1631	-9.5019     -13.4668    -11.1115	-26.1074 ;...
                    -20.7620	-7.8400     -11.1115	-9.1681 	-21.5412 ;...
                    -48.7823	-18.4208	-26.1074   	-21.5412	-50.6129 ];
           
start_cov = 0;
k_cov = (Effects*(max(s(1:10)/(SDEN_2*abyv)-start_cov,0)/(1.0-start_cov)));
A6_Cov = [0;k_cov(2:6);0;k_cov(2:6)]/R_e;
end

