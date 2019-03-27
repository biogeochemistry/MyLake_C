function [z_epi_obs] = z_epilimnion_obs(Tzt,zz)

Tzt = Tzt';
% Calculates the thickess of the epilimnion from the temperature
% observations, i.e. the last depth that exist in the epilimnion

z_epi_obs = -1*ones(1,size(Tzt,2));

epil_thres=1.0;

for k = 1:size(Tzt,2)
    Tz = Tzt(:,k);
    dS2dz = abs(Tz(1)-Tz);
    dk=find((dS2dz>epil_thres));
    %try
    if(sum(isnan(Tz))>0)
        z_epi_obs(k) = NaN;
    elseif(~isempty(dk))
        z_epi_obs(k) = zz(dk(1));
    end
    if(k<0)
        keyboard
    end
    %if(k==241);keyboard;end
end

z_epi_obs(z_epi_obs==-1) = max(zz);
end