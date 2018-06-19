function [CO2_new,DIC_new,Hplus_new] = DIC_diffusion2(DIC_old,Tz,Hplus_old,Fi_O2)
%Calculates new IC equilibrium after diffusion.

global ies80;

density = (polyval(ies80,max(0,Tz))+min(Tz,0))*0.001; %kg/l

DIC_old = 0.001.*DIC_old; %mg/l

%Carbon acid solubility and dissociation constants (Millero, 1995)
[K1, K2, Kw] = dissconst(Tz); % mol/kg ; mol/kg; mol^2/kg^2; 

Hplus_old = Hplus_old./density; %mol/kg

CO2mfrac_old = Hplus_old.*Hplus_old./((Hplus_old.*Hplus_old+Hplus_old.*K1+K1.*K2)); %mol CO2 / mol DIC
HCO3mfrac_old = Hplus_old.*K1./((Hplus_old.*Hplus_old+Hplus_old.*K1+K1.*K2)); %mol HCO3 / mol DIC
CO3mfrac_old = K1.*K2./((Hplus_old.*Hplus_old+Hplus_old.*K1+K1.*K2)); %mol CO3 / mol DIC

M_DIC_old = CO2mfrac_old*44.01+HCO3mfrac_old*61.01+CO3mfrac_old*60.01; % DIC molar mass g/mol

% Old molal DIC concentration 
DICm_old = DIC_old./M_DIC_old; %DIC concentration in mmol DIC/l
DIC_old_m = 0.001*DICm_old./density; %DIC concentration in mol/kg

OHm_old = Kw./Hplus_old.*density; %(mol/kg)^2 * (mol/kg)^(-1) * kg/l = mol/l
ALm_old = (HCO3mfrac_old + 2*CO3mfrac_old).*DICm_old+OHm_old*1000-(Hplus_old.*density)*1000; %mmol/l = mol/m3
% Old molal alkalinity
%AL_old_m = ALm_old./(1000*density); % mol/kg

%Mola_r_ DIC diffuses
DICm_new = Fi_O2 \ DICm_old; % mmol/l

%Mola_r_ alkalinity diffuses
ALm_new = Fi_O2 \ ALm_old; % mmol/l

DIC_new_m = 0.001*DICm_new./density; %mol/kg
AL_new_m = ALm_new./(1000*density); %mol/kg

F = DIC_new_m; % mol/kg
AL = AL_new_m; % mol/kg

f4 = -1*ones(length(F),1);
f3 = -K1-AL;
f2 = Kw+K1.*(F-K2-AL);
f1 = K1.*(Kw+K2.*(2*F-AL));
f0 = Kw.*K1.*K2;

Hplus_new = NaN*ones(length(Hplus_old),1);

for k = 1:length(Hplus_old)
    dummy = roots([f4(k),f3(k),f2(k),f1(k),f0(k)]);
    ratios = abs(dummy/Hplus_old(k)-1);
    Hplus_new(k) = dummy(ratios==min(ratios)); %mol /kg
end

CO2mfrac_new = Hplus_new.*Hplus_new./((Hplus_new.*Hplus_new+Hplus_new.*K1+K1.*K2));
HCO3mfrac_new = Hplus_new.*K1./((Hplus_new.*Hplus_new+Hplus_new.*K1+K1.*K2));
CO3mfrac_new = K1.*K2./((Hplus_new.*Hplus_new+Hplus_new.*K1+K1.*K2));

CO2m_new = CO2mfrac_new.*DICm_new; %mmol/l (= mol/m3)
HCO3m_new = HCO3mfrac_new.*DICm_new; % (mol HCO3 / mol DIC) * (mol DIC / m3) = mmol HCO3 / l
CO3m_new = CO3mfrac_new.*DICm_new; %mmol/l

M_DIC_new = CO2mfrac_new*44.01+HCO3mfrac_new*61.01+CO3mfrac_new*60.01; % new DIC molar mass g/mol

DIC_new = 1000*DICm_new.*M_DIC_new; %mg/g * mol/m3 * g/mol = mg/m3
CO2_new = 1000*44.01*CO2m_new; % mg/g * g/mol * mol/m3 = mg/m3

%Uusi alkaliniteetti
OHm_new = (Kw./Hplus_new).*density; %mol/l
ALm_uusi = (HCO3m_new + 2*CO3m_new)+OHm_new*1000-(Hplus_new.*density)*1000; %mmol/l
Hplus_new = Hplus_new.*density; %mol / l
%keyboard
end