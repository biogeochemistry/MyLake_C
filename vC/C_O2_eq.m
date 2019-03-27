function [O2_eq,K0_O2] = C_O2_eq(Pa,T0)
%Dissolved oxygen equilibrium concentration

%Inputs
%Pa          Air pressure (mbar = hPa)  
%T0          Water surface temperature

%Outputs
%O2_eq      DO equilibrium concentration (mg/m^3)

V_m = 22.414; % (L O2 / mol O2) O2 molar volume at STP (273.15 K, 1 atm)

Pa = 0.98692e-3*Pa; %atm
T0 = max(T0,0);
O2_ppm = 210000; %O2 concentration in air (mumol O2 / mol air)

if(T0 < 0)
    T0 = 0;
end    

%O2 solubility constant
lnbeeta = -58.3877+85.8079*(100./(T0+273.15))+23.8439*log((T0+273.15)./100);
beeta = exp(lnbeeta); %Bunsen solubility coefficient (L O2 / (L water * atm))
K0_O2 = beeta./V_m; %(mol O2 / (L water * atm))

%O2 equilibrium concentration (g/mol * mol/(L*atm) * mumol/mol * atm) = mug/l = mg/m^3 
O2_eq = 32*K0_O2.*O2_ppm.*Pa;

end