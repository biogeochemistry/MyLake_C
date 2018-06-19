
% Module for calculating river inflow and its effect on profiles of IC system

% INPUTS:
%   vertical arrays:    z is model grid, Vz is layer volume, Tz layer property (e.g. temperature);
%   scalars:            lvlD is the grid depth level *above* which the inflow settles (m), 
%                       Iflw is the inflow volume (m3/day), and 
%                       T_Iflw is the property of the inflow, and dz the grid stepsize (layer thickness)
% OUTPUTS: Cz is the new property profile after inflow

function [CO2_new,DIC_new,Hplus_new]=IOflow_DIC(~, zz, Vz, Tz, DIC_old, Hplus_old, lvlD, Iflw, Iflw_T, Iflw_DIC, Iflw_Hplus)

global ies80;
density = (polyval(ies80,max(0,Tz))+min(Tz,0))*0.001; %kg/l
density_iflw = (polyval(ies80,max(0,Iflw_T))+min(Iflw_T,0))*0.001; %kg/l

dum=find(zz>=lvlD); %same as lvlD but in grid level numbers
lvlG=dum(1);

%Inflow properties
%Tarvitaan inflown moolimäärä - ja alkaliniteetti, mitä sille tarkkaan
%ottaen sitten tapahtuukaan...

Iflw_DIC = 0.001.*Iflw_DIC; %mg/l

%Carbon acid solubility and dissociation constants (Millero, 1995)
[K1, K2, Kw] = dissconst(Iflw_T); % mol/kg ; mol/kg; mol^2/kg^2; 

Iflw_Hplus = Iflw_Hplus./density_iflw; %mol/kg

CO2mfrac_iflw = Iflw_Hplus.*Iflw_Hplus./((Iflw_Hplus.*Iflw_Hplus+Iflw_Hplus.*K1+K1.*K2)); %mol CO2 / mol DIC
HCO3mfrac_iflw = Iflw_Hplus.*K1./((Iflw_Hplus.*Iflw_Hplus+Iflw_Hplus.*K1+K1.*K2)); %mol HCO3 / mol DIC
CO3mfrac_iflw = K1.*K2./((Iflw_Hplus.*Iflw_Hplus+Iflw_Hplus.*K1+K1.*K2)); %mol CO3 / mol DIC

M_DIC_iflw = CO2mfrac_iflw*44.01+HCO3mfrac_iflw*61.01+CO3mfrac_iflw*60.01; % DIC molar mass g/mol

% Old molar DIC concentration in the unstable layer
DICm_iflw = Iflw_DIC./M_DIC_iflw; %DIC concentration in mmol DIC/l
DIC_iflw_m = 0.001*DICm_iflw./density_iflw; %DIC concentration in mol/kg

OHm_iflw = Kw./Iflw_Hplus.*density_iflw; %(mol/kg)^2 * (mol/kg)^(-1) * kg/l = mol/l
ALm_iflw = (HCO3mfrac_iflw + 2*CO3mfrac_iflw).*DICm_iflw+OHm_iflw*1000-(Iflw_Hplus.*density_iflw)*1000; %mmol/l = mol/m3
% Old molal alkalinity in the unstable layer
AL_iflw_m = ALm_iflw./(1000*density_iflw); % mol/kg

%Vastaavasti vesikerroksen DIC-moolimäärä ja alkaliniteetti

DIC_old = 0.001.*DIC_old; %mg/l

%Carbon acid solubility and dissociation constants (Millero, 1995)
[K1, K2, Kw] = dissconst(Tz); % mol/kg ; mol/kg; mol^2/kg^2; 

Hplus_old = Hplus_old./density; %mol/kg

CO2mfrac_old = Hplus_old.*Hplus_old./((Hplus_old.*Hplus_old+Hplus_old.*K1+K1.*K2)); %mol CO2 / mol DIC
HCO3mfrac_old = Hplus_old.*K1./((Hplus_old.*Hplus_old+Hplus_old.*K1+K1.*K2)); %mol HCO3 / mol DIC
CO3mfrac_old = K1.*K2./((Hplus_old.*Hplus_old+Hplus_old.*K1+K1.*K2)); %mol CO3 / mol DIC

M_DIC_old = CO2mfrac_old*44.01+HCO3mfrac_old*61.01+CO3mfrac_old*60.01; % DIC molar mass g/mol

% Old molar DIC concentration in the unstable layer
DICm_old = DIC_old./M_DIC_old; %DIC concentration in mmol DIC/l
DIC_old_m = 0.001*DICm_old./density; %DIC concentration in mol/kg

OHm_old = Kw./Hplus_old.*density; %(mol/kg)^2 * (mol/kg)^(-1) * kg/l = mol/l
ALm_old = (HCO3mfrac_old + 2*CO3mfrac_old).*DICm_old+OHm_old*1000-(Hplus_old.*density)*1000; %mmol/l = mol/m3
% Old molal alkalinity in the unstable layer
AL_old_m = ALm_old./(1000*density); % mol/kg

%Koska tasojen lämpötilat eivät taustalla muutu, mmol/l = mol/m3 sekoittuu.

if (lvlD==0) %if inflow is lighter than surface water, mix it with first layer
    DICm_new = DICm_old;
    DICm_new(1)=(Vz(1)*DICm_old(1)+Iflw*DICm_iflw)/(Vz(1)+Iflw);
    ALm_new = ALm_old;
    ALm_new(1)=(Vz(1)*ALm_old(1)+Iflw*ALm_iflw)/(Vz(1)+Iflw);
    
else    %otherwise add the inflow to the appropriate
    %depth level and "lift" the water column above
    DICm_new = DICm_old;
    DICrev_pour=flipud([DICm_old(1:lvlG-1); DICm_iflw; 0]); %up-side down Tz to be poured in (plus an extra zero)
    ALm_new = ALm_old;
    ALrev_pour=flipud([ALm_old(1:lvlG-1); ALm_iflw; 0]);
    Vzrev_pour=flipud([Vz(1:lvlG-1); Iflw; 0]);  %up-side down Vz to be poured in (plus an extra zero)
    Vzrev_fill=flipud(Vz(1:lvlG-1));  %up-side down Vz to be filled in
    
    for i=1:length(Vzrev_fill)
        inxA=find(cumsum(Vzrev_pour)>Vzrev_fill(i));
        inxB=inxA(1); %index of the last layer to be poured in (partly)
        ShakerV=[Vzrev_pour(1:inxB-1); Vzrev_fill(i)-sum(Vzrev_pour(1:inxB-1))];
        ShakerDIC=DICrev_pour(1:inxB);
        ShakerAL=ALrev_pour(1:inxB);
        ALm_new(lvlG-i)=sum(ShakerV.*ShakerAL)/sum(ShakerV); %new property after mixing
        DICm_new(lvlG-i)=sum(ShakerV.*ShakerDIC)/sum(ShakerV); %new property after mixing
        Vzrev_pour(1:inxB)=Vzrev_pour(1:inxB)-ShakerV; %subtract poured volumes from the "reserves"
    end
    
    
end %if

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
uusi_OHm = (Kw./Hplus_new)./density; %mol/l
alkalinitym_uusi = (HCO3m_new + 2*CO3m_new)+uusi_OHm*1000-(Hplus_new)*1000; %mmol/l
Hplus_new = Hplus_new.*density; %mol / l
%keyboard
end