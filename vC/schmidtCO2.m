function Sc = schmidtCO2(Tw)
%Sc = schmidt(Tw) computes the Schmidt number for CO2 as a
%function of temperature using Jähne et al. (1987)
%
%Inputs
%Tw     Water temperature (oC)
%
%Outputs
%Sc    Schmidt number for CO2(-)

T = Tw;

A = 1911.1; %Schmidt number polynomial fit coefficients
B = 118.11;
C = 3.4527;
D = 0.041320;

Sc = A-B*T+C*T.^2-D*T.^3;

end