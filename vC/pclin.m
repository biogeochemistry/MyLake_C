function [PCzt,TCzt,mlt] = pclin(Tzt,zz)

% Calculate pycnocline depth
dz = zz(2)-zz(1);

TCzt = NaN*ones(1,size(Tzt,2));
PCzt = NaN*ones(1,size(Tzt,2));
mlt = NaN*ones(1,size(Tzt,2));

pycno_thres=0.1;  %threshold density gradient value (kg m-3 m-1)
thermo_thres=1.0;
epil_thres=thermo_thres;
ies80 = [6.536332e-9,-1.120083e-6,1.001685e-4,-9.09529e-3,6.793952e-2,999.842594];
for k = 1:size(Tzt,2)
    Tz = Tzt(:,k);
    rho = polyval(ies80,max(0,Tz(:))) + min(Tz(:),0);
    dRdz = [NaN; abs(diff(rho))];
    di=find((dRdz<(pycno_thres*dz)) | isnan(dRdz));
    dRdz(di)=0;
    PCz = sum(zz .* dRdz) ./ sum(dRdz);
    PCzt (k) = PCz;
    
    dSdz = [NaN; abs(diff(Tz))];
    dj=find((dSdz<(thermo_thres*dz)) | isnan(dSdz));
    dSdz(dj)=0;
    TCz = sum(zz .* dSdz) ./ sum(dSdz);
    TCzt(k) = TCz;
    
    dS2dz = abs(Tz(1)-Tz);
    dk=find((dS2dz>epil_thres));
    %try
        if(~isempty(dk))
            mlt(k) = zz(dk(1));
        end
    if(k<0)
        keyboard
    end
    %if(k==159);keyboard;end
end
%Muutos DIC-version määritelmään
PCzt(isnan(PCzt)) = max(zz);
TCzt(isnan(TCzt)) = max(zz);
mlt(isnan(mlt)) = max(zz);
end