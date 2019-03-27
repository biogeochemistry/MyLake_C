function [CO2_new,surfaceflux,CO2_eq,K0,ma] = CO2flux(k_CO2,CO2_old,Pa,Ts,dz,currdate,Az,Vz,ppm_base)

%Inputs
%k_CO2       Gas exchange coefficient (m/s)
%CO2_old     Old CO2 concentration in the surface layer (mg/m^3) 
%Pa          Air pressure (mbar = hPa)  
%Ts          Water surface temperature
%dz          Surface layer depth
%currdate    Date number
%Az          Water layer area (m^3)
%Vz          Water layer volume (m^3)
%ppm_base    CO2 atmospheric mixing ratio (mumol CO2 / mol air)

%Outputs
%CO2_new     New CO2 concentration (mg/m^3)
%surfaceflux CO2 air-water surface flux mg/(m^2*d)
%CO2_eq      Carbon dioxide equilibrium concentration (mg/m^3)
%K0          Carbon dioxide solubility coefficient (mol/(l*atm))
%CO2_ppm     Carbon dioxide concentration in air (mumol CO2 / mol air)

CO2_new = CO2_old;

CO2_air_ppm = (ppm_base+4*sin(((currdate+0.310625)/365.2425*2*pi))); 

%CO2 equilibrium concentration (mg/m^3)

[CO2_eq, K0] = C_CO2_eq(CO2_air_ppm,Pa,Ts);

alpha = 1; %Chemical enhancement factor

surfaceflux = 86400*alpha*k_CO2*(CO2_old(1)-CO2_eq); %(s/d * m/s * mg/m^3 = mg/(m^2*d))

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
            mass_wc = mass_wc+(CO2_old(3)-CO2_eq)*Vz(3);
            levels = [1 2 3];
            if(mass_wc < mass)
                mass_wc = mass_wc+(CO2_old(4)-CO2_eq)*Vz(4);
                levels = [1 2 3 4];
                if(mass_wc < mass)
                    mass_wc = mass_wc+(CO2_old(5)-CO2_eq)*Vz(5);
                    levels = [1 2 3 4 5];
                    if(mass_wc < mass)
                        mass_wc = mass_wc+(CO2_old(6)-CO2_eq)*Vz(6);
                        levels = [1 2 3 4 5 6];
                        if(mass_wc < mass)
                            mass_wc = mass_wc+(CO2_old(7)-CO2_eq)*Vz(7);
                            levels = [1 2 3 4 5 6 7];
                            if(mass_wc < mass)
                                %mass_wc = mass_wc+(CO2_old(8)-CO2_eq)*Vz(8);
                                levels = [1 2 3 4 5 6 7 8];
                                disp(['CO2-vuo menee kahdeksanteen kerrokseen päivänä ', num2str(currdate)]);
                            end
                        end
                    end
                end
            end
        end
    end
    
    %Riittävien kerrosten hiilidioksidin kokonaismassa miinus pois menevä.
    mass_final= sum(Vz(levels).*CO2_old(levels))-mass;
    CO2_new(levels) = mass_final./sum(Vz(levels)); %(mg/m^3); time step = 1 d)
    ma = length(levels);
    if(CO2_new(1)<CO2_eq);disp(['CO2-vuo imee liikaa päivänä ', num2str(currdate)]);end
    %if(CO2_old(1)/CO2_eq <1.3);keyboard;end
else
    CO2_new(1) = CO2_old(1)-surfaceflux/dz; %(mg/m^3); time step = 1 d)
    ma = -1;
end
%CO2_new = max(0,CO2_old-surfaceflux/dz);
%if(CO2_old-surfaceflux/dz < 0)
%    surfaceflux = (CO2_new - CO2_old)*dz;
%    keyboard
%end
%if(currdate>735784);keyboard;end
if(CO2_new(1)<0);keyboard;end
end
