function alpha = thermexp(Tw)
%alpha = thermexp(Tw) computes the thermal expansion coefficient of water as a
%function of temperature using ies80 as water density
%
%Inputs
%Tw     Water temperature (oC)
%
%Outputs
%thermexp    Water thermal expansion coefficient (1/K)

%Water density polynomial
global ies80;

T = max(0,Tw)+min(Tw,0);

rho = polyval(ies80,T); %kg/m^3

drho = 6.793952e-2-1.819058e-2*T+3.005055e-4*T.*T-4.480332e-6*T.^3+3.268166e-8*T.^4; %kg/(m^3*T)

alpha = -drho./rho;