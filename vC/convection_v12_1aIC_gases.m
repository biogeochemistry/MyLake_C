% === MyLake model, version 1.2, 15.03.05 ===
% by Tom Andersen & Tuomo Saloranta, NIVA 2004

% VERSION 1.2.1a, based on convection_v12 (with three modified lines of code, marked with NEW!!!)

% Convection module
% Code checked by TSA, 07.03.05
% Last modified by TSA, 17.07.07
% Modified by PK 30.12.2010 (DIC) & 14.02.2011 (O2) & 03.11.2015 (POC)

function [Tz,Cz,Sz,Pz,Chlz,PPz,DOPz,DOCz1,DOCz2,DOCz3,DICz,O2z,POCz1,POCz2,Hplusz,CO2z] = ...
    convection_v12_1aIC_gases(Tz_in,Cz_in,Sz_in,Pz_in,Chlz_in,PPz_in,DOPz_in,DOCz1_in,DOCz2_in,DOCz3_in,DICz_in,O2z_in,POCz1_in,POCz2_in,Hplusz_in,CO2z_in,Tprof_prev,Vz,Cw,f_par,lambdaz_wtot_avg,zz,swa_b0,tracer_switch,springautumn,i)

% Inputs (with extension "_in") and Outputs:
%       Tz   : Temperature profile
%       Cz   : Tracer profile
%       Sz   : Suspended inorg. matter profile
%       Pz   : Dissolved inorg. P profile
%       Chlz : Chlorophyll a profile
%       PPz  : Phosphorus bound to inorganic particles profile
%       DOPz  : Dissolved organic phosphorus profile
%       DOCz  : Particulate inorganic phosphorus profile
%       DICz  : Dissolved inorganic carbon profile (PK)
%       O2z  : Dissolved oxygen profile (PK)
%       POCz  : Particulate organic carbon profile (PK)
%       Hplusz : Hydronium ion concentration

% Inputs:
%       Tprof_prev   : Temperature profile from previous timestep
%       etc.

% These variables are still global and not transferred by functions
global ies80;
%T_dummy = Tz_in;
Trhomax=3.98; %temperature of maximum water density (deg C)
Nz=length(zz); %total number of layers in the water column
dz=zz(2)-zz(1); %model grid step

%pres = sw_pres(zz+dz/2,61.5);

% Convective mixing adjustment
% Mix successive layers until stable density profile
rho = polyval(ies80,max(0,Tz_in(:)))+min(Tz_in(:),0);	% Density (kg/m3)
%rho = sw_dens(Salz_in,max(0,Tz_in(:)),0) + min(Tz_in(:),0);
d_rho=[diff(rho); 1]; %d_rho = how much a layer is lighter than layer below; last cell in "d_rho" is always positive (sediment)
if(i==-1);keyboard;end
while any(d_rho < 0),
    blnUnstb_layers=(d_rho <= 0); %1=layer is heavier or equal than layer below, 0=layer is lighter than layer below
    A_Unstb=find(diff([0; blnUnstb_layers])==1); %layer index(es) where unstable/neutral column(s) start(s)
    B_Unstb=find(diff([0; blnUnstb_layers])==-1)-1;%layer index(es) where unstable/neutral column(s) end(s)

    for n = 1:length(A_Unstb)
        j = [A_Unstb(n):B_Unstb(n)+1];
        T_old = Tz_in(j);
        Tmix = sum(Tz_in(j) .* Vz(j)) / sum(Vz(j));
        Tz_in(j) = Tmix * ones(size(Tz_in(j)));

        if (tracer_switch==1)
            Cmix = sum(Cz_in(j) .* Vz(j)) / sum(Vz(j));
            Cz_in(j) = Cmix * ones(size(Cz_in(j)));
        end

        DIC_old = DICz_in(j);
        Hplus_old = Hplusz_in(j);
        T_new = Tz_in(j);
        
        [CO2z_in(j),DICz_in(j),Hplusz_in(j)] = DIC_convection2(DIC_old,T_old,T_new,Hplus_old,Vz(j));
        
        Smix = sum(Sz_in(j) .* Vz(j)) / sum(Vz(j));
        Sz_in(j) = Smix * ones(size(Sz_in(j)));

        Pmix = sum(Pz_in(j) .* Vz(j)) / sum(Vz(j));
        Pz_in(j) = Pmix * ones(size(Pz_in(j)));

        Chlmix = sum(Chlz_in(j) .* Vz(j)) / sum(Vz(j));
        Chlz_in(j) = Chlmix * ones(size(Chlz_in(j)));

        PPmix = sum(PPz_in(j) .* Vz(j)) / sum(Vz(j));
        PPz_in(j) = PPmix * ones(size(PPz_in(j)));

        DOPmix = sum(DOPz_in(j) .* Vz(j)) / sum(Vz(j));
        DOPz_in(j) = DOPmix * ones(size(DOPz_in(j)));

        DOC1mix = sum(DOCz1_in(j) .* Vz(j)) / sum(Vz(j));
        DOCz1_in(j) = DOC1mix * ones(size(DOCz1_in(j)));
        
        DOC2mix = sum(DOCz2_in(j) .* Vz(j)) / sum(Vz(j));
        DOCz2_in(j) = DOC2mix * ones(size(DOCz2_in(j)));
        
        DOC3mix = sum(DOCz3_in(j) .* Vz(j)) / sum(Vz(j));
        DOCz3_in(j) = DOC3mix * ones(size(DOCz3_in(j)));
        
        O2mix = sum(O2z_in(j) .* Vz(j)) / sum(Vz(j));
        O2z_in(j) = O2mix * ones(size(O2z_in(j)));
        
        POC1mix = sum(POCz1_in(j) .* Vz(j)) / sum(Vz(j));
        POCz1_in(j) = POC1mix * ones(size(POCz1_in(j)));
        
        POC2mix = sum(POCz2_in(j) .* Vz(j)) / sum(Vz(j));
        POCz2_in(j) = POC2mix * ones(size(POCz2_in(j)));
        
    end
    if(i==-1);keyboard;end
    rho = polyval(ies80,max(0,Tz_in(:))) + min(Tz_in(:),0);
    %rho = sw_dens(Salz_in,max(0,Tz_in(:)),0) + min(Tz_in(:),0);
    d_rho=[diff(rho); 1];
end;

if (springautumn==1)
    % Spring/autumn turnover
    % don't allow temperature jumps over temperature of maximum density
    Tz_interm = Tz_in;
    if( ((Tprof_prev(1)>Trhomax)&(Tz_in(1)<Trhomax))|((Tprof_prev(1)<Trhomax)&(Tz_in(1)>Trhomax)) ) %NEW!!! (Tprof, ">" instead of ">=")!

        jumpinx=find( ((Tprof_prev>Trhomax)&(Tz_in<Trhomax))|((Tprof_prev<Trhomax)&(Tz_in>Trhomax)) ); %NEW!!!!!
        if (sum(jumpinx==1)==0);disp('NOTE: Non-surface jumps over temperature of maximum density');end %NEW!!!!!

        intSign=sign(Trhomax-Tz_in(1)); %plus in autumn turnover, minus in spring turnover
        XE_turn=cumsum((Tz_in-Trhomax).*Vz*Cw*intSign); %always starts negative
        Dummy=find(XE_turn>0);
        if(isempty(Dummy)==1)
            Tz_in(:)=Trhomax;
            %O2mix = sum(O2z_in .* Vz) / sum(Vz);
            %O2z_in(:) = O2mix;
            if (intSign==1)
                Tz_in(1)=Tz_in(1)+intSign*XE_turn(end)/(Vz(1)*Cw); %put overshoot on top layer
            else
                Tz_in=Tz_in + (-diff( intSign*XE_turn(end) * (f_par * exp([0; -lambdaz_wtot_avg] .* [zz; zz(end)+dz]) + ...
                    (1-f_par) * exp([0; -swa_b0*ones(Nz,1)] .* [zz; zz(end)+dz])) ) ./(Vz*Cw));
                %distribute overshoot as shortwave energy
            end
        else
            Tz_in(1:Dummy(1)-1)=Trhomax;
            Tz_in(Dummy(1))=Trhomax+intSign*XE_turn(Dummy(1))/(Vz(Dummy(1))*Cw);
        end
        
      [CO2z_in,DICz_in,Hplusz_in] = DICsystem_new_gases(DICz_in,CO2z_in,Tz_interm,Tz_in,Hplusz_in,i)  ;
      
      %Voi olla tilanne, että jäätä vasten oleva kerros ei lämpene yli 4
      %asteen mutta seuraava lämpenee, jolloin konvektio estyy.
      
%     elseif(Tprof_prev(1)==0)%Jään alla tulee vain kerran springautumn-ehto
%         
%         if( ((Tprof_prev(2)>Trhomax)&&(Tz_in(2)<Trhomax))||((Tprof_prev(2)<Trhomax)&&(Tz_in(2)>Trhomax)) ) %NEW!!! (Tprof, ">" instead of ">=")!
%             keyboard
%             jumpinx=find( ((Tprof_prev(2:end)>Trhomax)&(Tz_in(2:end)<Trhomax))|((Tprof_prev(2:end)<Trhomax)&(Tz_in(2:end)>Trhomax)) ); %NEW!!!!!
%             if (sum(jumpinx==1)==0);disp('NOTE: Non-surface jumps over temperature of maximum density');end %NEW!!!!!
%             
%             intSign=sign(Trhomax-Tz_in(2)); %plus in autumn turnover, minus in spring turnover
%             XE_turn=cumsum((Tz_in(2:end)-Trhomax).*Vz(2:end)*Cw*intSign); %always starts negative
%             Dummy=find(XE_turn>0);
%             if(isempty(Dummy)==1)
%                 Tz_in((2:end))=Trhomax;
%                 if (intSign==1)
%                     Tz_in(2)=Tz_in(2)+intSign*XE_turn(end)/(Vz(2)*Cw); %put overshoot on top layer
%                 else
%                     Tz_in(2:end)=Tz_in(2:end) + (-diff( intSign*XE_turn(end) * (f_par * exp([0; -lambdaz_wtot_avg(2:end)] .* [zz(2:end); zz(end)+dz]) + ...
%                         (1-f_par) * exp([0; -swa_b0*ones(Nz-1,1)] .* [zz(2:end); zz(end)+dz])) ) ./(Vz(2:end)*Cw));
%                     %distribute overshoot as shortwave energy
%                 end
%             else
%                 Tz_in(2:Dummy(1)-1)=Trhomax;
%                 Tz_in(Dummy(1))=Trhomax+intSign*XE_turn(Dummy(1))/(Vz(Dummy(1))*Cw);
%             end
%             keyboard
%             [CO2z_in,DICz_in,Hplusz_in] = DICsystem_new(DICz_in,CO2z_in,Tz_interm,Tz_in,Hplusz_in,i);
%         end
    end
    
    
end %springautumn
if(i==-116);keyboard;end
Tz=Tz_in;
Cz=Cz_in;
Sz=Sz_in;
Pz=Pz_in;
PPz=PPz_in;
Chlz=Chlz_in;
DOPz=DOPz_in;
DOCz1=DOCz1_in;
DOCz2=DOCz2_in;
DOCz3=DOCz3_in;
DICz=DICz_in;
O2z=O2z_in;
POCz1=POCz1_in;
POCz2=POCz2_in;
Hplusz=Hplusz_in;
CO2z=CO2z_in;if(i==-1);keyboard;end
%keyboard
end