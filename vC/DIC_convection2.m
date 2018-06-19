function [CO2_new,DIC_new,Hplus_new] = DIC_convection2(DIC_old,T_old,T_new,Hplus_old,Vz)
%Calculates new pH due to the change in CO2 concentration.

global ies80;

density_old = (polyval(ies80,max(0,T_old))+min(T_old,0))*0.001; %kg/l
density_new = (polyval(ies80,max(0,T_new))+min(T_new,0))*0.001; %kg/l

DIC_old = 0.001.*DIC_old; %mg/l
%T_old = max(0,T_old)+273.15;
T_new = max(0,T_new)+273.15;

%Carbon acid solubility and dissociation constants (Millero, 1995)
% K1 = 290.9097-14554.21./T_old-45.0575*log(T_old); %~mol/kg% K1 = exp(K1); %HCO3; mol/kg
% K2 = 207.6548-11843.79./T_old-33.6485*log(T_old); %~mol/kg% K2 = exp(K2); %CO3; mol/kg
% Kw = 148.9802-13847.26./T_old-23.6521*log(T_old); %~mol/kg% Kw = exp(Kw); %H2O; (mol/kg)^2
[K1, K2, Kw] = dissconst(T_old); % mol/kg ; mol/kg; mol^2/kg^2; 

Hplus_old = Hplus_old./density_old; %mol/kg

CO2mfrac_old = Hplus_old.*Hplus_old./((Hplus_old.*Hplus_old+Hplus_old.*K1+K1.*K2)); %mol CO2 / mol DIC
HCO3mfrac_old = Hplus_old.*K1./((Hplus_old.*Hplus_old+Hplus_old.*K1+K1.*K2)); %mol HCO3 / mol DIC
CO3mfrac_old = K1.*K2./((Hplus_old.*Hplus_old+Hplus_old.*K1+K1.*K2)); %mol CO3 / mol DIC

M_DIC_old = CO2mfrac_old*44.01+HCO3mfrac_old*61.01+CO3mfrac_old*60.01; % DIC molar mass g/mol

% Old molar DIC concentration in the unstable layer
DICm_old = DIC_old./M_DIC_old; %DIC concentration in mmol DIC/l
DIC_old_m = 0.001*DICm_old./density_old; %DIC concentration in mol/kg

OHm_old = Kw./Hplus_old.*density_old; %(mol/kg)^2 * (mol/kg)^(-1) * kg/l = mol/l
alkalinitym_old = (HCO3mfrac_old + 2*CO3mfrac_old).*DICm_old+OHm_old*1000-(Hplus_old.*density_old)*1000; %mmol/l = mol/m3
% Old molal alkalinity in the unstable layer
AL_old_m = alkalinitym_old./(1000*density_old); % mol/kg

%Water volume to water mass
Mz_old = 1000*Vz.*density_old; % l/m3* m3 * kg/l = kg
DIC_moles = sum(DIC_old_m .* Mz_old); %mol
AL_moles = sum(AL_old_m .* Mz_old); %mol
Mz_new =  1000*Vz.*density_new; %kg

DICmix = DIC_moles / sum(Mz_new); % mol/kg
% New molal DIC concentration in the mixed zone
DIC_new_m = DICmix * ones(size(DIC_old_m)); % mol/kg

ALmix = AL_moles / sum(Mz_new); % mol/kg
% New molal alkalinity in the mixed zone
AL_new_m = ALmix * ones(size(AL_old_m)); % mol/kg

DICm_new = 1000*DIC_new_m.*density_new; %mmol/l

F = DIC_new_m; % mol/kg
AL = AL_new_m;

K1 = 290.9097-14554.21./T_new-45.0575*log(T_new); %~mol/kg
K1 = exp(K1); %HCO3; mol/kg
K2 = 207.6548-11843.79./T_new-33.6485*log(T_new); %~mol/kg
K2 = exp(K2); %CO3; mol/kg
Kw = 148.9802-13847.26./T_new-23.6521*log(T_new); %~mol/kg
Kw = exp(Kw); %H2O; (mol/kg)^2

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
HCO3z_new = HCO3m_new*61.01;
CO3z_new = CO3m_new*60.01;

%Uusi alkaliniteetti
uusi_OHm = (Kw./Hplus_new)./density_new; %mol/l
alkalinitym_uusi = (HCO3m_new + 2*CO3m_new)+uusi_OHm*1000-(Hplus_new)*1000; %mmol/l
Hplus_new = Hplus_new.*density_new; %mol / l
%keyboard
end