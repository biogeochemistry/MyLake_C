function [CO2z_new,DICz_new,Hplus_new] = DICsystem_new(DICz,CO2z_in,Tz_old,Tz_new,Hplus,i)
%Calculates new pH due to the change in CO2 concentration.
%Inputs
%Hplus       Old 10^(-pH)
%DICz        Old DIC profile (mg/m^3)
%CO2z_in     CO2 profile after production/consumption processes (mg/m^3)
%Tz          water temperature profile (deg C)

%Outputs
%CO2z_new    carbon dioxide (mg/m^3)
%DICz_new    dissolved inorganic carbon (mg/m^3)
%Hplus_new   10^(-pH)

%Water density polynomial 
%if(i==511);keyboard;end
global ies80;
%ies80 = [6.536332e-9,-1.120083e-6,1.001685e-4,-9.09529e-3,6.793952e-2,999.842594];
density_old = (polyval(ies80,max(0,Tz_old))+min(Tz_old,0))*0.001; %kg/l
density_new = (polyval(ies80,max(0,Tz_new))+min(Tz_new,0))*0.001; %kg/l
% Note: in equations of density it is assumed that every supercooled degree lowers density by 
% 1 kg m-3 due to frazil ice formation (probably no practical meaning, but included for "safety")

DICz = 0.001.*DICz; %mg/l
Tz_old = max(0,Tz_old)+273.15;
Tz_new = max(0,Tz_new)+273.15;

%Carbon acid solubility and dissociation constants (Millero, 1995)

K1 = 290.9097-14554.21./Tz_old-45.0575*log(Tz_old); %~mol/kg
K1 = exp(K1); %HCO3; mol/kg
K2 = 207.6548-11843.79./Tz_old-33.6485*log(Tz_old); %~mol/kg
K2 = exp(K2); %CO3; mol/kg
Kw = 148.9802-13847.26./Tz_old-23.6521*log(Tz_old); %~mol/kg
Kw = exp(Kw); %H2O; (mol/kg)^2

%Hydrogen ion concentration

Hplus = Hplus./density_old; %mol/kg

%Old fractions of carbonates in a particular form
%http://www.chem.usu.edu/~sbialkow/Classes/3650/Carbonate/Carbonic%20Acid.html

CO2mfrac = Hplus.*Hplus./((Hplus.*Hplus+Hplus.*K1+K1.*K2)); %mol CO2 / mol DIC
HCO3mfrac = Hplus.*K1./((Hplus.*Hplus+Hplus.*K1+K1.*K2)); %mol HCO3 / mol DIC
CO3mfrac = K1.*K2./((Hplus.*Hplus+Hplus.*K1+K1.*K2)); %mol CO3 / mol DIC

M_DIC = CO2mfrac*44.01+HCO3mfrac*61.01+CO3mfrac*60.01; % DIC molar mass g/mol

DICmz = DICz./M_DIC; %DIC concentration in mmol DIC/l

%vanha_DIC_mooleissa = DICmz;

CO2mz_old = CO2mfrac.*DICmz; % mmol/l

%CO2z_old = 44.01.*CO2mfrac.*DICmz; % (g CO2 / mol CO2) * (mol CO2 / mol DIC) * (mol DIC / m3) = mg CO2 / l
%HCO3z = 61.01.*HCO3mfrac.*DICmz; % (g HCO3 / mol HCO3) * (mol HCO3 / mol DIC) * (mol DIC / m3) = mg HCO3 / l
%CO3z = 60.01.*CO3mfrac.*DICmz; % (g CO3 / mol CO3) * (mol CO3 / mol DIC) * (mol DIC / m3) = mg CO3 / l

%HCO3mz = HCO3mfrac.*DICmz; % (mol HCO3 / mol DIC) * (mol DIC / m3) = mmol HCO3 / l
%CO3mz = CO3mfrac.*DICmz; %(mol CO3 / mol DIC) * (mol DIC / m3) = mmol CO3 / l
OHm = Kw./Hplus.*density_old; %(mol/kg)^2 * (mol/kg)^(-1) * kg/l = mol/l
alkalinitym_vanha = (HCO3mfrac + 2*CO3mfrac).*DICmz+OHm*1000-(Hplus.*density_old)*1000; %mmol/l = mol/m3
AL = alkalinitym_vanha./(1000*density_new); % mol/kg

%Hiilidioksidin moolim‰‰r‰ muuttuu, mik‰ muuttaa myˆs DICin moolim‰‰r‰‰ 
%(tasapainossa kolme osaa muuttuvat toisikseen, joten tasapainohaussa moolim‰‰r‰ ei muutu)

CO2z_in = 0.001*CO2z_in; %mg CO2 / l

CO2mz = (1/44.01)*CO2z_in; % (mmol CO2 / mg CO2) * (mg CO2 / l) = mmol CO2 / l

CO2z_m_intermediate = 0.001*CO2mz./density_new; % mol/kg

CO2z_m_old = 0.001*CO2mz_old./density_new; % mol/kg

DIC_old_m = 0.001*DICmz./density_new; % mol/kg

DIC_new_m = DIC_old_m + (CO2z_m_intermediate-CO2z_m_old); % mol/kg

DICmz_new = 1000*DIC_new_m.*density_new; %mmol/l

F = DIC_new_m; % mol/kg

K1 = 290.9097-14554.21./Tz_new-45.0575*log(Tz_new); %~mol/kg
K1 = exp(K1); %HCO3; mol/kg
K2 = 207.6548-11843.79./Tz_new-33.6485*log(Tz_new); %~mol/kg
K2 = exp(K2); %CO3; mol/kg
Kw = 148.9802-13847.26./Tz_new-23.6521*log(Tz_new); %~mol/kg
Kw = exp(Kw); %H2O; (mol/kg)^2

f4 = -1*ones(length(F),1);
f3 = -K1-AL;
f2 = Kw+K1.*(F-K2-AL);
f1 = K1.*(Kw+K2.*(2*F-AL));
f0 = Kw.*K1.*K2;

Hplus_new = NaN*ones(length(Hplus),1);
%if(i==511);keyboard;end
for k = 1:length(Hplus)
    dummy = roots([f4(k),f3(k),f2(k),f1(k),f0(k)]);
    dummy(dummy<0) = NaN;
    ratios = abs(dummy/Hplus(k)-1);
    try
        Hplus_new(k) = dummy(ratios==min(ratios)); %mol /kg
    catch
        keyboard
    end
end

CO2mfrac_new = Hplus_new.*Hplus_new./((Hplus_new.*Hplus_new+Hplus_new.*K1+K1.*K2));
HCO3mfrac_new = Hplus_new.*K1./((Hplus_new.*Hplus_new+Hplus_new.*K1+K1.*K2));
CO3mfrac_new = K1.*K2./((Hplus_new.*Hplus_new+Hplus_new.*K1+K1.*K2));

CO2mz_new = CO2mfrac_new.*DICmz_new; %mmol/l
HCO3mz_new = HCO3mfrac_new.*DICmz_new; % (mol HCO3 / mol DIC) * (mol DIC / m3) = mmol HCO3 / l
CO3mz_new = CO3mfrac_new.*DICmz_new; %mmol/l

M_DIC_new = CO2mfrac_new*44.01+HCO3mfrac_new*61.01+CO3mfrac_new*60.01; % new DIC molar mass g/mol

DICz_new = 1000*DICmz_new.*M_DIC_new; %mg/m3
CO2z_new = 1000*44.01*CO2mz_new; % mg/g * g/mol * mol/m3 = mg/m3
%HCO3z_new = HCO3mz_new*61.01;
%CO3z_new = CO3mz_new*60.01;

%Uusi alkaliniteetti
uusi_OHm = (Kw./Hplus_new)./density_new; %mol/l
alkalinitym_uusi = (HCO3mz_new + 2*CO3mz_new)+uusi_OHm*1000-(Hplus_new.*density_new)*1000; %mmol/l
Hplus_new = Hplus_new.*density_new; %mol / l
%if(~isreal(Hplus_new));keyboard;end
%if(i==181);keyboard;end
end

% ekat = NaN*ones(length(Hplus),4);
% tokat = NaN*ones(length(Hplus),4);
% kolkit = NaN*ones(length(Hplus),1);
% 
% for k = 1:length(Hplus)
%     ekat(k,:) = roots([f4(k),f3(k),f2(k),f1(k),f0(k)]);
%     ekat(k,(ekat(k,:)<0)) = NaN;
%     tokat(k,:) = abs(ekat(k,:)/Hplus(k)-1);
%     kolkit(k) = ekat(k,tokat(k,:)==min(tokat(k,:))); %mol /kg
%     %keyboard
% end


