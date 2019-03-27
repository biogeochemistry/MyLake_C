function visc = viscwater(Tw)
%vis = viscwater(Tw) computes the kinematic viscosity of water as a
%function of temperature following White (2011)
%
%Inputs
%Tw     Temperature (oC)
%
%Outputs
%visc    Kinematic viscosity (m^2/s)

%Water density polynomial
global ies80;

visc_dyn_0 = 1.788e-3; %kg/(m*s)

density = polyval(ies80, max(0,Tw)+min(Tw,0)); %kg/m^3

z = 273./(Tw+273.15);
x = -1.704-5.306.*z+7.003.*z.*z;
visc_dyn = visc_dyn_0*exp(x);

visc = visc_dyn./density;

