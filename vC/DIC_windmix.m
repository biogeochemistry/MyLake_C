function [CO2_new,DIC_new,Hplus_new] = DIC_windmix(DIC_old,T_old,T_new,Hplus_old,Vz,zb,KP_ratio)
%Calculates IC wind mixing.

global ies80;

density_old = (polyval(ies80,max(0,T_old))+min(T_old,0))*0.001; %kg/l
density_new = (polyval(ies80,max(0,T_new))+min(T_new,0))*0.001; %kg/l

DIC_old = 0.001.*DIC_old; %mg/l

%Carbon acid solubility and dissociation constants (Millero, 1995)
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

Mz_old = 1000*Vz.*density_old; % l/m3* m3 * kg/l = kg
Mz_new = 1000*Vz.*density_new; %kg

if(KP_ratio>=1) %if all water is mixed
    DIC_moles = sum(DIC_old_m .* Mz_old); %mol
    AL_moles = sum(AL_old_m .* Mz_old); %mol
    DICmix = DIC_moles / sum(Mz_new); %mol/kg
    % New molal DIC concentration in the mixed zone
    DIC_new_m=DICmix * ones(size(DIC_old_m)); % mol/kg
    ALmix = AL_moles / sum(Mz_new); % mol/kg
    % New molal alkalinity in the mixed zone
    AL_new_m = ALmix * ones(size(AL_old_m)); % mol/kg
else %if the layers 1:zb are fully mixed and the layer zb+1 partially mixed
    DIC_new_m = DIC_old_m; %mol/kg
    AL_new_m = AL_old_m; %mol/kg
    DICmix = sum( [Mz_old(1:zb); KP_ratio*Mz_old(zb+1)].*DIC_old_m(1:zb+1) ) / sum([Mz_new(1:zb); KP_ratio*Mz_new(zb+1)]);
    DIC_new_m(1:zb) = DICmix; %mol/kg
    DIC_new_m(zb+1) = (KP_ratio*DICmix*Mz_new(zb+1) + (1-KP_ratio)*DIC_old_m(zb+1)*Mz_old(zb+1)) ./ Mz_new(zb+1);
    ALmix = sum( [Mz_old(1:zb); KP_ratio*Mz_old(zb+1)].*AL_old_m(1:zb+1) ) / sum([Mz_new(1:zb); KP_ratio*Mz_new(zb+1)]);
    AL_new_m(1:zb) = ALmix; %mol/kg
    AL_new_m(zb+1)= (KP_ratio*ALmix*Mz_new(zb+1) + (1-KP_ratio)*AL_old_m(zb+1)*Mz_old(zb+1)) ./ Mz_new(zb+1);
    
end

DICm_new = 1000*DIC_new_m.*density_new; %mmol/l

F = DIC_new_m; % mol/kg
AL = AL_new_m;

[K1, K2, Kw] = dissconst(T_new);

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
alkalinitym_uusi = (HCO3m_new + 2*CO3m_new)+uusi_OHm*1000-(Hplus_new.*density_new)*1000; %mmol/l
Hplus_new = Hplus_new.*density_new; %mol / l
%keyboard
end