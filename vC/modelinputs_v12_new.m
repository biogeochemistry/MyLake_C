% === MyLake model, version 1.2, 15.03.05 ===
% by Tom Andersen & Tuomo Saloranta, NIVA 2005
%
% Module for reading input data and parameters
% Code checked by TSA, xx.03.2005
% Last modified by TSA, 15.08.2006 (Az replaced by In_Az 10.03.06; Possibility to have NaN in Global rad. series, 15.08.06)          
% Modified by PK 30.12.2010 (DIC) & 11.2.2011 (oxygen) & 3.11.2015 (POC) & 27.7.2016 (F_NLOM + pH)

function [In_Z,In_Az,tt,In_Tz,In_Cz,In_Sz,In_TPz,In_DOPz,In_Chlz,In_DOCz,In_DICz,In_O2z,In_POCz,In_TPz_sed,In_FNLOM,In_FIM,In_pH,Ice0,OC_frac0,Wt,Inflw,...
         Phys_par,Phys_par_range,Phys_par_names,Bio_par,Bio_par_range,Bio_par_names] ...
            = modelinputs_v12_new(M_start,M_stop,init_filename,init_sheet,...
            input_filename,input_sheet,param_filename,param_sheet,dt)

% Inputs:
%       M_start : Model start date [year, month, day]
%       M_stop : Model stop date [year, month, day]
%       + Input filenames and sheetnames
%		dt		: Time step (= 1 day)
% Outputs:
%		tt		: Solution time domain (day)
%       In_Z    : Depths read from initial profiles file (m)
%       In_Az   : Areas read from initial profiles file (m2)
%       In_Tz   : Initial temperature profile read from initial profiles file (deg C)
%       In_Cz   : Initial tracer profile read from initial profiles file (-)
%       In_Sz   : Initial sedimenting tracer (or suspended inorganic matter) profile read from initial profiles file (kg m-3)
%       In_TPz  : Initial total P profile read from initial profiles file (mg m-3)
%       In_DOPz  : Initial dissolved organic P profile read from initial profiles file (mg m-3)
%       In_Chlz : Initial chlorophyll a profile read from initial profiles file (mg m-3)
%       In_DOCz  : Initial DOC profile read from initial profiles file (mg m-3)
%       In_DICz  : Initial DIC profile read from initial profiles file (mg m-3) (PK)
%       In_O2z  : Initial oxygen profile read from initial profiles file (mg m-3) (PK)
%       In_POCz  : Initial POC profile read from initial profiles file (mg m-3) (PK)
%       In_TPz_sed  : Initial total P profile in the sediment compartments read from initial profiles file (mg m-3)
%       In_NLOM      : Initial profile of volume fraction of nonliving organic matter in the sediment solids (dry weight basis)
%       In_FIM      : Initial profile of volume fraction of inorganic matter in the sediment solids (dry weight basis)
%       In_pH      : Initial pH profile read from initial profiles file (-)
%       Ice0            : Initial conditions, ice and snow thicknesses (m) (Ice, Snow)
%		Wt		        : Weather data
%       Inflow          : Inflow data
%       Phys_par        : Main 23 parameters that are more or less fixed
%       Phys_par_range  : Minimum and maximum values for Phys_par (23 * 2)
%       Phys_par_names  : Names for Phys_par
%       Bio_par         : Main 15 parameters  that are more or less site specific
%       Bio_par_range   : Minimum and maximum values for Bio_par (15 * 2)
%       Bio_par_names   : Names for Bio_par

global ies80;

% == Read model parameter file
[ParaMx,StrMx]=xlsread(param_filename,param_sheet);

%Main physical parameters (dz, Kz_ak, etc...)	
Phys_par_names=StrMx(3:25,1);
Phys_par=ParaMx(3:25,2);
Phys_par_range=ParaMx(3:25,3:4);

%Main biological parameters (Y_cp, m_twty, g_twty, etc...)	
Bio_par_names=StrMx(26:65,1);
Bio_par=ParaMx(26:65,2);
Bio_par_range=ParaMx(26:65,3:4);

  % Vertical settling velocities
  U = Bio_par(8:9); %for sedimenting velocities
  if any(U<0)
   error('Given settling velocity must be positive')  
  end
  
% == Read morphometric and initial profile file

 [InitMx,StrMx]=xlsread(init_filename,init_sheet);
 In_Z=InitMx(3:end,1);
 In_Az=InitMx(3:end,2);
 In_Tz=InitMx(3:end,3);
 In_Cz=InitMx(3:end,4);
 In_Sz=InitMx(3:end,5);
 In_TPz=InitMx(3:end,6);
 In_DOPz=InitMx(3:end,7);
 In_Chlz=InitMx(3:end,8);
 In_DOCz=InitMx(3:end,9);
 In_DICz=InitMx(3:end,10);
 In_O2z=InitMx(3:end,11);
 In_POCz=InitMx(3:end,12);
 In_TPz_sed=InitMx(3:end,13);
 In_FNLOM=InitMx(3:end,14);
 In_FIM=InitMx(3:end,15);
 In_pH=InitMx(3:end,16);
 
 Ice0=InitMx(3,17:18);
 OC_frac0=InitMx(3,19:21);
 
   tt = [datenum(M_start):dt:datenum(M_stop)]';		% Solution time domain

% == Read input forcing data file

[InputMx,StrMx]=xlsread(input_filename,input_sheet);

In_Date=InputMx(3:end,1:3);
In_Met=InputMx(3:end,4:10);
In_Inflow=InputMx(3:end,11:22);

tmet=datenum(In_Date);

if(tmet(1) > datenum(M_start))
    error('Simulation start date is too early')
elseif(tmet(end) < datenum(M_start))
    error('Simulation start date is too late')
elseif(tmet(end) < datenum(M_stop))
    error('Simulation end date is too late')
end

dum=100*((tmet(end)-tmet(1)+1)-length(tmet))/(tmet(end)-tmet(1)+1);
disp('Percent missing dates in meteorology and inflow data: ');
disp([num2str(dum) ' %']);

dum=100*sum(isnan(In_Met))./length(tmet);
disp('Percent missing values in meteorology data (values correspond to columns 4-10 in input file): ');
disp([num2str(dum) ' %']);

dum=100*sum(isnan(In_Inflow))./length(tmet);
disp('Percent missing values in inflow data (values correspond to columns 11-21 in input file): ');
disp([num2str(dum) ' %']);
disp(' ')

clear Wt
for i=1:7 %Interpolate over missing values and dates
    nonnans = find(isnan(In_Met(:,i))==0);
    if(isempty(nonnans)) % if the whole column is NaNs then preserve it
        Wt(:,i) = NaN*ones(length(tt(:)),1);   
    else
        repaired = interp1(nonnans,In_Met(nonnans,i),[1:length(In_Met(:,i))]);
        Wt(:,i) = interp1(tmet, repaired, tt(:));
    end
end
% Wt(:,1)  Global radiation (MJ/(m^2 day))
% Wt(:,2)  Cloud cover (-)
% Wt(:,3)  Air temperature (deg. C, at 2 m height)
% Wt(:,4)  Relative humidity (%, at 2 m height)
% Wt(:,5)  Air pressure (mbar)
% Wt(:,6)  Wind speed (m/s at 10 m height)
% Wt(:,7)  Precipitation (mm/day)

clear Inflw
for i=1:12 %Interpolate over missing values and dates
    nonnans = find(isnan(In_Inflow(:,i))==0);
        if(isempty(nonnans)) % if the whole column is NaNs then preserve it
         Inflw(:,i) = NaN*ones(length(tt(:)),1);   
        else
         repaired = interp1(nonnans,In_Inflow(nonnans,i),[1:length(In_Inflow(:,i))]);
         Inflw(:,i) = interp1(tmet, repaired, tt(:));
        end
 end
% Inflw(:,1) Inflow volume (m3 day-1)
% Inflw(:,2) Inflow temperature (deg C)
% Inflw(:,3) Inflow tracer concentration (-)
% Inflw(:,4) Inflow sedimenting tracer (or suspended inorganic matter) concentration (kg m-3)
% Inflw(:,5) Inflow total phosphorus (TP) concentration (mg m-3)
% Inflw(:,6) Inflow dissolved organic phosphorus (DOP) concentration (mg m-3)
% Inflw(:,7) Inflow chlorophyll a concentration (mg m-3)
% Inflw(:,8) Inflow DOC concentration (mg m-3)
% Inflw(:,9) Inflow DIC concentration (mg m-3) (PK)
% Inflw(:,10) Inflow O2 concentration (mg m-3) (PK)
% Inflw(:,11) Inflow POC concentration (mg m-3) (PK)
% Inflw(:,12) Inflow H+ concentration (mg m-3) (PK)

% International Equation of State 1980
% 5-order polynomial for density as function of temperature
ies80 = [6.536332e-9,-1.120083e-6,1.001685e-4,-9.09529e-3,6.793952e-2,999.842594];


% Default turbulence and wind shelter parameterization (Hondzo and Stefan, 1993; Ellis et al., 1991)
if(isnan(Phys_par(2)))
  Phys_par(2) = 0.00706*(In_Az(1)/1e6)^0.56; % default diffusion coeff. parameterisation
end

if(isnan(Phys_par(3)))
  Phys_par(3) = 8.98e-4;		%default value for diffusion coeff. in ice-covered water   
end

if(isnan(Phys_par(4)))
  Phys_par(4) = 7e-5;			% default minimum allowed stability frequency, N2 > N0 <=> Kz < Kmax (1/s2)    		
end

if(isnan(Phys_par(5)))
  Phys_par(5) =  1-exp(-0.3*In_Az(1)/1e6);			% default wind sheltering parameterisation		
end
