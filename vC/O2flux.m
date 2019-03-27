function [O2_new,surfaceflux,O2_eq,K0_O2,ma] = O2flux(k_O2,O2_old,Pa,Ts,currdate,Az,Vz)

%Inputs
%k_O2        Gas exchange coefficient (m/s)
%O2_old      Old O2 concentration (mg/m^3) 
%Pa          Air pressure (mbar = hPa)  
%Ts          Water surface temperature
%dz          Surface layer depth
%Az          Water layer area (m^3)
%Vz          Water layer volume (m^3)

%Outputs
%O2_new      New O2 concentration (mg/m^3)
%surfaceflux O2 air-water surface flux mg/(m^2*d)
%O2_eq       O2 equilibrium concentration (mg/m^3)
%K0_O2       O2 solubility coefficient (mol/(l*atm))


O2_new = O2_old;

%O2 equilibrium concentration (mg/m^3)
[O2_eq, K0_O2] = C_O2_eq(Pa,Ts);

surfaceflux = 86400*k_O2*(O2_old(1)-O2_eq); %(s/d * m/s * mg/m^3 = mg/(m^2*d))

if(surfaceflux > 0)
    %Lasketaan vaihtuva kokonaismassa
    mass = surfaceflux*Az(1);
    %Ylimmän kerroksen hapen kokonaismassa
    mass_wc = (O2_old(1)-O2_eq)*Vz(1);
    levels = 1;
    %Jos tämä ei riitä,
    if(mass_wc < mass)
        %otetaan mukaan toinen kerros.
        mass_wc = mass_wc+(O2_old(2)-O2_eq)*Vz(2);
        levels = [1 2];
        if(mass_wc < mass)
            mass_wc = mass_wc+(O2_old(3)-O2_eq)*Vz(3);
            levels = [1 2 3];
            %disp(['Menee kolmanteen kerrokseen päivänä ', num2str(date)])
            if(mass_wc < mass)
                mass_wc = mass_wc+(O2_old(4)-O2_eq)*Vz(4);
                levels = [1 2 3 4];
                %disp(['Menee neljänteen kerrokseen päivänä ', num2str(date)]);
                if(mass_wc < mass)
                    levels = [1 2 3 4 5];
                    %disp(['Menee viidenteen kerrokseen päivänä ', num2str(date)]);
                end
            end
        end
    end
    
    %Riittävien kerrosten hapen kokonaismassa miinus pois menevä.
    mass_final= sum(Vz(levels).*O2_old(levels))-mass;
    O2_new(levels) = mass_final./sum(Vz(levels)); %(mg/m^3); time step = 1 d)
    ma = length(levels);
    %if(O2_new(1)<O2_eq);disp(['Vuo imee liikaa päivänä ', num2str(date)]);keyboard;end
else
    %O2_new(1) = O2_old(1)-surfaceflux.*Az(1)./Vz(1); %(mg/m^3); time step = 1 d)
    
    %Lasketaan vaihtuva kokonaismassa
    mass = -surfaceflux*Az(1);
    %Ylimmän kerroksen hapen kokonaismassa
    mass_wc = -(O2_old(1)-O2_eq)*Vz(1);
    levels = 1;
    %Jos tämä ei riitä,
    if(mass_wc < mass)
        %otetaan mukaan toinen kerros.
        mass_wc = mass_wc-(O2_old(2)-O2_eq)*Vz(2);
        levels = [1 2];
        if(mass_wc < mass)
            mass_wc = mass_wc-(O2_old(3)-O2_eq)*Vz(3);
            levels = [1 2 3];
            if(mass_wc < mass)
                mass_wc = mass_wc-(O2_old(4)-O2_eq)*Vz(4);
                levels = [1 2 3 4];
                if(mass_wc < mass)
                    mass_wc = mass_wc-(O2_old(5)-O2_eq)*Vz(5);
                    levels = [1 2 3 4 5];
                    if(mass_wc < mass)
                        mass_wc = mass_wc-(O2_old(6)-O2_eq)*Vz(6);
                        levels = [1 2 3 4 5 6];
                        if(mass_wc < mass)
                            mass_wc = mass_wc-(O2_old(7)-O2_eq)*Vz(7);
                            levels = [1 2 3 4 5 6 7];
                            if(mass_wc < mass)
                                %mass_wc = mass_wc-(O2_old(8)-O2_eq)*Vz(8);
                                levels = [1 2 3 4 5 6 7 8];
                                disp(['Happivuo menee kahdeksanteen kerrokseen päivänä ', num2str(currdate)]);
                                
                            end
                        end
                    end
                end
            end
        end
    end
    %if(length(levels)>3);keyboard;end
    mass_final= sum(Vz(levels).*O2_old(levels))+mass;
    O2_new(levels) = mass_final./sum(Vz(levels)); %(mg/m^3); time step = 1 d)
    ma = -length(levels);
end
%if(O2_new(1)<0);keyboard;end
end
