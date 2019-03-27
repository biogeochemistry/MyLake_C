function Fb = buoyancyflux(Heff, Tw)
%Fb = buoyancyflux(Heff, Tw) computes buoyancy flux using effective heat
%flux (Imberger 1985) and mixing layer temperature
%
%Inputs
%Heff   Effective heat flux (W/m^2)
%Tw     Water temperature (oC)
%
%Outputs
%Fb     Buoyancy flux (m^2/s^2)

%Water density polynomial
global ies80;

g = 9.81; %(m/s^2)
cp = 4.18e3; % Specific heat capacity of water (MyLake default) (J/(kg*K))

rho_w = polyval(ies80, max(0,Tw)+min(Tw,0)); %kg/m^3
alpha = thermexp(Tw);
try
Fb = g*alpha.*Heff./(cp*rho_w);
catch
    keyboard
end
end