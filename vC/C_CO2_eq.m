function [CO2_eq,K0] = C_CO2_eq(CO2_air,Pa,T0)
%Carbon dioxide equilibrium concentration

%Inputs
%CO2_air     Air CO2 concentration (ppm)
%Pa          Air pressure (mbar = hPa)  
%T0          Water surface temperature

%Outputs
%CO2_eq      CO2 equilibrium concentration (mg/m^3)


global ies80;

density = (polyval(ies80,max(0,T0))+min(T0,0))*0.001; %kg/l

Pa = 0.98692e-3*Pa; %atm
T0 = max(T0,0);

%CO2 solubility constant

K0 = -60.2409+93.4517.*(100./(T0+273.15))+23.3585*log((T0+273.15)/100); %~mol/(kg*atm)
K0 = exp(K0).*density; %CO_2; mol/(kg*atm)*kg/l = mol/(l*atm)

%CO2 equilibrium concentration (g/mol*mol/(l*atm)*(mumol/mol)*atm = mug/l = mg/m^3)

CO2_eq = 44.01.*K0.*CO2_air.*Pa;

end
