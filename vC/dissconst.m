function [K1,K2,Kw] = dissconst(Tz)
%Calculates inorganic carbon dissosiation constants.

Tz = max(0,Tz)+273.15;

%Carbon acid solubility and dissociation constants (Millero, 1995)

K1 = 290.9097-14554.21./Tz-45.0575*log(Tz); %~mol/kg
K1 = exp(K1); %HCO3; mol/kg
K2 = 207.6548-11843.79./Tz-33.6485*log(Tz); %~mol/kg
K2 = exp(K2); %CO3; mol/kg
Kw = 148.9802-13847.26./Tz-23.6521*log(Tz); %~mol/kg
Kw = exp(Kw); %H2O; (mol/kg)^2
end