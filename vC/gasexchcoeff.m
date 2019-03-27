function [k,dissip] = gasexchcoeff(Heff, ustar_a, Tw, Ws, zML, Ta, Rh, Pa, gas, k_equation)
%gec = gasexchcoeff(Heff, Tw) computes buoyancy flux using effective heat
%
%
%Inputs
%Heff   Effective heat flux (W/m^2)
%Tw     Water temperature (oC)
%k_equation
%
%Outputs
%k          Gas exchange coefficient (m/s)

%Water density polynomial
global ies80;
if(sum(isnan(Rh))>0)
    Rh(isnan(Rh)) = 75;
    disp('NaN- suht. kosteuksia korvattu 75 %:lla')
end

if(sum(isnan(Pa))>0)
    Pa(isnan(Pa)) = 1000;
    disp('NaN-ilmanpaineita korvattu 1000 Pa:lla')
end

switch gas
    case 'CO2'
        Sc = schmidtCO2(Tw);
    case 'O2'
        Sc = schmidtO2(Tw);
end


switch k_equation
    case 'Heiskanen'
        
        C1 = 0.00015;
        C2 = 0.07;
        
        Fb = buoyancyflux(Heff, Tw);
        
        u_ref = C1*Ws;
        wstar_w = zeros(size(Fb));
        for m = 1:length(wstar_w)
            if(Fb(m)>0)
                wstar_w(m) = 0;
            else
                wstar_w(m) = (-Fb(m).*zML(m)).^(1/3);
            end
        end
        k = sqrt((u_ref).^2+(C2*wstar_w).^2).*Sc.^(-0.5);
        dissip = [];
        %if(abs(Ws-5.3962)<0.01);keyboard;end
    case 'SR_MacIntyre'
        
        c1 = 0.50;
        c2 = 0.77;
        c3 = 0.3;
        kappa = 0.4;
        z = 0.15;
        
        rho_water = polyval(ies80, Tw);
        rho_air = air_dens(Ta,Rh,Pa);
        ustar_w = ustar_a.*sqrt(rho_air./rho_water);
        
        Fb = buoyancyflux(Heff, Tw);
        nu_w = viscwater(Tw);
        
        k = zeros(size(Fb));
        dissip = zeros(size(Fb));
        for m = 1:length(k)
            dissip(m) = c3*ustar_w(m).^3/(kappa*z)-c2*Fb(m);
            if(dissip(m)<0)
                k(m) = c1*((c3*ustar_w(m).^3/(kappa*z)-0*Fb(m)).*nu_w(m)).^0.25.*Sc(m).^(-0.5); %m/s
                disp('SR:n dissipaatio menee pakkaselle solussa');keyboard
                disp(num2str(m))
            else
                k(m) = c1*((c3*ustar_w(m).^3/(kappa*z)-c2*Fb(m)).*nu_w(m)).^0.25.*Sc(m).^(-0.5); %m/s
            end
        end
            %keyboard
    case 'SR_Tedford'
        
        c1 = 0.5;
        c2 = 0.56;
        c3 = 0.77;
        c4 = 0.6;
        kappa = 0.4;
        z = 0.15;
        
        rho_water = polyval(ies80, Tw);
        rho_air = air_dens(Ta,Rh,Pa);
        ustar_w = ustar_a.*sqrt(rho_air./rho_water);
        
        Fb = buoyancyflux(Heff, Tw);
        nu_w = viscwater(Tw);
        
        dissip = zeros(size(Fb));
        for m = 1:length(Fb)
            if(Fb(m)<0)
                dissip(m) = c2*ustar_w(m).^3/(kappa*z)-c3*Fb(m);
            else
                dissip(m) = c4*ustar_w(m).^3/(kappa*z);
            end
        end
        k = c1*(dissip.*nu_w).^(0.25).*Sc.^(-0.5);
        %keyboard
        %if(abs(Pa-984.0902)<0.001);keyboard;end
    case 'CC'
        
        z0 = exp((1.22*log(1.5)-log(10))/0.22); %roughness length (m) fixed so that U(10) = 1.22*U(1.5) (~ 0.00027)
        z_meas = 1.5; %wind speed measurement height (m)
        Ufactor = (log(10)-log(z0))/(log(z_meas)-log(z0));
        U10 = Ufactor*Ws;
        
        k = (2.07+0.215*U10.^1.7).*(Sc/600).^(-0.5); %cm/h
        k = k/(3600*100); %m/s
        dissip = [];
        
    case 'U_MacIntyre'
                
        z0 = exp((1.22*log(1.5)-log(10))/0.22); %roughness length (m) fixed so that U(10) = 1.22*U(1.5) (~ 0.00027)
        z_meas = 1.5; %wind speed measurement height (m)
        Ufactor = (log(10)-log(z0))/(log(z_meas)-log(z0));
        U10 = Ufactor*Ws;
        
        Fb = buoyancyflux(Heff, Tw);
        
        k = zeros(size(U10));
        for m = 1:length(k)
            if(isnan(Fb(m)))
                k(m) = NaN;
            elseif(Fb(m)<0)
                k(m) = 2.04*U10(m)+2.0; %cm/h
            else
                k(m) = 1.74*U10(m)-0.15;
            end
        end
        k = k.*(Sc/600).^(-0.5)/(3600*100); %m/s
        dissip = [];
        %keyboard
    otherwise
        k = [];
        dissip = [];
        disp('Wrong gas exchange formula!')
        keyboard
end


end