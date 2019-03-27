function [Sc] = schmidtO2(Tw)
%Sc = schmidt(Tw) computes the Schmidt number for O2 as a
%function of temperature using the citation in Wanninkhof (1992)
%
%Inputs
%Tw     Water temperature (oC)
%
%Outputs
%Sc    Schmidt number for O2(-)

T = Tw;

A = 1800.6; %Schmidt number polynomial fit coefficients
B = 120.10;
C = 3.7818;
D = 0.047608;

Sc = A-B*T+C*T.^2-D*T.^3;
end