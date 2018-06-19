function [CO2_new,surfaceflux,CO2_eq,K0,CO2_ppm] = carbondioxideflux_new(CO2_old,Ws,Pa,T0,dz,date,Az,Vz,ppm_base)
%Carbon dioxide surface flux
%Yritetään äkkiä "hätäratkaisuna" uutta: arvot jäävät ihan lian pieniksi!
%Tai ainakin ensiksi CO":n meno pakkaselle. Ilmeisesti niin kävi
%Valkea-Kotisessakin, mutta koska siellä ei ole pH-laskua, negatiivinen DIC
%ja CO2 tasoittuivat wind mixingissä. Nyt CO2:sta lasketaan heti Hplus.
% -> No niin. Nyt pitäisi olla vähän paremmin vuo ulos, niin ettei mene
% pakkaselle.

%Inputs
%CO2_old     Old CO2 concentration in the surface layer (mg/m^3) 
%Ws          Wind speed (m/s)
%Pa          Air pressure (mbar = hPa)  
%T0          Water surface temperature
%dz          Surface layer depth
%date        Date number
%ppm_base    CO2 concentration in air (mumol CO2 / mol air)

%Outputs
%CO2_new     New CO2 concentration (mg/m^3)
%surfaceflux CO2 air-water surface flux mg/(m^2*d)
%K0          Carbon dioxide solubility coefficient (mol/(l*atm))
%CO2_ppm     Carbon dioxide concentration in air (mumol CO2 / mol air)

%Water density polynomial 

global ies80;

CO2_new = CO2_old;

density = (polyval(ies80,max(0,T0))+min(T0,0))*0.001; %kg/l
% Note: in equations of density it is assumed that every supercooled degree lowers density by 
% 1 kg m-3 due to frazil ice formation (probably no practical meaning, but included for "safety")

Pa = 0.98692e-3*Pa; %atm
T0 = max(T0,0);

%CO2 solubility constant

K0 = -60.2409+93.4517.*(100./(T0+273.15))+23.3585*log((T0+273.15)/100); %~mol/(kg*atm)
K0 = exp(K0)*density; %CO_2; mol/(kg*atm)*kg/l = mol/(l*atm)

%CO2 concentration in air (mumol CO2 / mol air)
%ppm_base = 380;
%default (kalibraatioarvo!) 380;
%1981-2010 362; 2071-2100 RCP2.6 429; 2071-2100 RCP4.5 533; 2071-2100 RCP8.5 807; 

CO2_ppm = (ppm_base+4*sin(((date+0.310625)/365.2425*2*pi))); 

%CO2 equilibrium concentration (g/mol*mol/(l*atm)*(mumol/mol)*atm = mug/l = mg/m^3)

CO2_eq = 44.01*K0*CO2_ppm*Pa;

alpha = 1; %Chemical enhancement factor
A = 1911.1; %Schmidt number polynomial fit coefficients
B = 118.11;
C = 3.4527;
D = 0.041320;

k_600 = 2.07+0.215*Ws^1.7; %transfer velocity for Schmidt number 600 (cm/h)
schmidt = A-B*T0+C*T0^2-D*T0^3; %Schmidt number for CO2 (-)
k_CO2 = k_600*(schmidt/600)^(-0.6667); %transfer velocity for CO2 (cm/h)

surfaceflux = 0.24*alpha*k_CO2*(CO2_old(1)-CO2_eq); %(m/cm * h/d * cm/h *mg/m^3 = mg/(m^2*d))
if(surfaceflux > 0)
    %Lasketaan vaihtuva kokonaismassa
    mass = surfaceflux*Az(1);
    %Ylimmän kerroksen hiilidioksidin kokonaismassa
    mass_wc = (CO2_old(1)-CO2_eq)*Vz(1);
    levels = 1;
    %Jos tämä ei riitä,
    if(mass_wc < mass)
        %otetaan mukaan toinen kerros.
        mass_wc = mass_wc+(CO2_old(2)-CO2_eq)*Vz(2);
        levels = [1 2];
        if(mass_wc < mass)
            %mass_wc = mass_wc+(CO2_old(3)-CO2_eq)*Vz(3);
            levels = [1 2 3];
        end
    end
    
    %Riittävien kerrosten hiilidioksidin kokonaismassa miinus pois menevä.
    mass_final= sum(Vz(levels).*CO2_old(levels))-mass;
    CO2_new(levels) = mass_final./sum(Vz(levels)); %(mg/m^3); time step = 1 d)
else
    CO2_new(1) = CO2_old(1)-surfaceflux/dz; %(mg/m^3); time step = 1 d)
end
%CO2_new = max(0,CO2_old-surfaceflux/dz);
%if(CO2_old-surfaceflux/dz < 0)
%    surfaceflux = (CO2_new - CO2_old)*dz;
%    keyboard
%end

end
