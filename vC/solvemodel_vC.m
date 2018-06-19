% === MyLake model, version 1.2, 15.03.05 ===
% by Tom Andersen & Tuomo Saloranta, NIVA 2005
%
% VERSION 1.2.1 (two phytoplankton groups are included; variable Cz denotes
% this second group now. Frazil ice included + some small bug-fixes and code rearrangements. Using convection_v12_1a.m code)
%
% Main module
% Code checked by TSA, xx.xx.200x
% Last modified by TSA, 21.08.2007

% Modified to include Fokema-module by Kai Rasmus. 16.5.2007
% Uusimman Fokema-version mukainen 30.12.2010
% Lisätty matriisit DOCzt1,DOCzt2,DOCzt3,Daily_BB1t,Daily_BB2t,Daily_BB3t,Daily_PBt
% jotka eivät nyt kuitenkaan lähde ulos funktiosta. -> Muokkaus 080915: lähtevät.

% Lisätty 29.12.2010 DIC-muuttuja, joka saa inflown ja sekoittuu
% konvektiolla ja diffuusiolla DOCin tapaan.

%Lisätään 10.2.2011 O2-muuttuja DICin mallin mukaan.

%Lisätään POC-muuttuja
%Ja siirretään xOC-dynamiikka samaan P-Chl-dynamiikan kanssa, siis niin että
%P-Chl-diffuusio tapahtuu samalla kuin xOC-diffuusio.

%Muokataan 27.7.2016 initejä ja parametrejä tiedostoista luettavaksi. Aiempi versio
%jää entiselleen. Tähän ei outputiin tule (vielä) muutoksia, vain
%funktion modelinputs_new.m kautta tulevaan kamaan.

%Tähän korjattu CO2-DIC-Hplus-tasapainoyhtälöt ja diffuusiot.

function [zz,Az,Vz,tt,Qst,Kzt,Tzt,Czt,Szt,Pzt,Chlzt,PPzt,DOPzt,DOCzt,DOCzt1,DOCzt2,DOCzt3,DOCtfrac,...
    Daily_BB1t,Daily_BB2t,Daily_BB3t,Daily_PBt,DICzt,CO2zt,O2zt,O2_sat_relt,O2_sat_abst,POCzt,BODzt,Qzt_sed,lambdazt,...
    P3zt_sed,P3zt_sed_sc,His,DoF,DoM,MixStat,Wt,surfaceflux,O2fluxt,CO2_eqt,K0t,O2_eqt,K0_O2t,...
    CO2_ppmt,dO2Chlt,POC1tfrac,dO2SODt,dO2DOCt,pHt,testi3t,T_sedt] = ...
    solvemodel_vC(M_start,M_stop,Initfile,Initsheet,Inputfile,Inputsheet,Parafile,Parasheet,Qlambda,varargin)

warning off MATLAB:fzero:UndeterminedSyntax %suppressing a warning message

% Inputs (to function)
%       M_start     : Model start date [year, month, day]
%       M_stop      : Model stop date [year, month, day]
%               + Input filenames and sheetnames

% Inputs (received from input module):
%		tt		: Solution time domain (day)
%       In_Z    : Depths read from initial profiles file (m)
%       In_Az   : Areas read from initial profiles file (m2)
%       In_Tz   : Initial temperature profile read from initial profiles file (deg C)
%       In_Cz   : Initial chlorophyll (group 2) profile read from initial profiles file (-)
%       In_Sz   : Initial sedimenting tracer (or suspended inorganic matter) profile read from initial profiles file (kg m-3)
%       In_TPz  : Initial total P profile read from initial profiles file (incl. DOP & Chla & Cz) (mg m-3)
%       In_DOPz  : Initial dissolved organic P profile read from initial profiles file (mg m-3)
%       In_Chlz : Initial chlorophyll (group 1) profile read from initial profiles file (mg m-3)
%       In_DOCz  : Initial DOC profile read from initial profiles file (mg m-3)
%       In_DICz  : Initial DIC profile read from initial profiles file (mg m-3) (PK)
%       In_O2z  : Initial oxygen profile read from initial profiles file (mg m-3) (PK)
%       In_POCz  : Initial POC profile read from initial profiles file (mg m-3) (PK)
%       In_TPz_sed  : Initial total P profile in the sediment compartments read from initial profiles file (mg m-3)
%       In_NLOM      : Initial profile of volume fraction of nonliving organic matter in the sediment solids (dry weight basis) (PK)
%       In_FIM      : Initial profile of volume fraction of inorganic matter in the sediment solids (dry weight basis)
%       In_pH      : Initial pH profile read from initial profiles file (-) (PK)
%       Ice0            : Initial conditions, ice and snow thicknesses (m) (Ice, Snow)
%		Wt		        : Weather data
%       Inflow          : Inflow data
%       Phys_par        : Main 23 parameters that are more or less fixed
%       Phys_par_range  : Minimum and maximum values for Phys_par (23 * 2)
%       Phys_par_names  : Names for Phys_par
%       Bio_par         : Main 23 parameters that are more or less site specific
%       Bio_par_range   : Minimum and maximum values for Bio_par (23 * 2)
%       Bio_par_names   : Names for Bio_par

% Outputs (other than Inputs from input module):
%		Qst : Estimated surface heat fluxes ([sw, lw, sl] * tt) (W m-2)
%		Kzt	: Predicted vertical diffusion coefficient (tt * zz) (m2 d-1)
%		Tzt	: Predicted temperature profile (tt * zz) (deg C)
%		Czt	: Predicted chlorophyll (group 2) profile (tt * zz) (-)
%		Szt	: Predicted passive sedimenting tracer (or suspended inorganic matter) profile (tt * zz) (kg m-3)=(g L-1)
%		Pzt	: Predicted dissolved inorganic phosphorus profile (tt * zz) (mg m-3)
%		Chlzt	    : Predicted chlorophyll (group 1) profile (tt * zz) (mg m-3)
%		PPzt	    : Predicted particulate inorganic phosphorus profile (tt * zz) (mg m-3)
%		DOPzt	    : Predicted dissolved organic phosphorus profile (tt * zz) (mg m-3)
%		DOCzt	    : Predicted dissolved organic carbon (DOC) profile (tt * zz) (mg m-3)
%		DICzt	    : Predicted dissolved inorganic carbon (DIC) profile (tt * zz) (mg m-3) (PK)
%		CO2zt	    : Predicted dissolved carbon dioxide profile (tt * zz) (mg m-3) (PK)
%		O2zt	    : Predicted dissolved oxygen profile (tt * zz) (mg m-3) (PK)
%       POCzt       : Predicted dissolved organic carbon (DOC) profile (tt * zz) (mg m-3)
%       O2_sat_rel  : Predicted relative oxygen saturation (PK)
%       O2_sat_abs  : Predicted absolute oxygen saturation (PK)
%		Qz_sed      : Predicted  sediment-water heat flux (tt * zz) (W m-2, normalised to lake surface area)
%       lambdazt    : Predicted average total light attenuation coefficient down to depth z (tt * zz) (m-1)
%       P3zt_sed    : Predicted P conc. in sediment for P (mg m-3), PP(mg kg-1 dry w.) and Chl (mg kg-1 dry w.) (tt * zz * 3)
%       P3zt_sed_sc : Predicted P source from sediment for P, PP and Chl (mg m-3 day-1) (tt * zz * 3)
%       His         : Ice information matrix ([Hi Hs Hsi Tice Tair rho_snow IceIndicator] * tt)
%       DoF, DoM    : Days of freezing and melting (model timestep number)
%       MixStat     : Temporary variables used in model testing, see code (N * tt)

% Fokema outputs
%       CDOMzt      : Coloured dissolved organic matter absorption m-1
%                   : (tt * zz)

% These variables are still global and not transferred by functions
global ies80;

tic
disp(['Running MyLake C from ' datestr(datenum(M_start)) ' to ' datestr(datenum(M_stop)) ' ...']);

% ===Switches===
snow_compaction_switch=1;       %snow compaction: 0=no, 1=yes
river_inflow_switch=1;          %river inflow: 0=no, 1=yes
sediment_heatflux_switch=1;     %heatflux from sediments: 0=no, 1=yes
selfshading_switch=1;           %light attenuation by chlorophyll a: 0=no, 1=yes
tracer_switch=1;                %simulate tracers:  0=no, 1=yes
DOC_attenuation_switch = 1;     %light attenuation by DOC: 0:no, 1=yes
density_switch = 0; %must be 0
% ==============

dt=1.0; %model time step = 1 day (DO NOT CHANGE!)

if (nargin>9) %if optional command line parameter input is used
    disp('Bypassing input files...Running with input data & parameters given on command line');
    [In_Z,In_Az,tt,In_Tz,In_Cz,In_Sz,In_TPz,In_DOPz,In_Chlz,In_DOCz,In_DICz,In_O2z,In_POCz,In_TPz_sed,In_FNLOM,In_FIM,In_pH,Ice0,OC_frac0,Wt,Inflw,...
        Phys_par,Phys_par_range,Phys_par_names,Bio_par,Bio_par_range,Bio_par_names]...
        = deal(varargin{:});
else
    %Read input data
    [In_Z,In_Az,tt,In_Tz,In_Cz,In_Sz,In_TPz,In_DOPz,In_Chlz,In_DOCz,In_DICz,In_O2z,In_POCz,In_TPz_sed,In_FNLOM,In_FIM,In_pH,Ice0,OC_frac0,Wt,Inflw,...
        Phys_par,Phys_par_range,Phys_par_names,Bio_par,Bio_par_range,Bio_par_names]...
        = modelinputs_v12_new(M_start,M_stop,Initfile,Initsheet,Inputfile,Inputsheet,Parafile,Parasheet,dt);
end

load albedot1.mat; %load albedot1 table, in order to save execution time

% Unpack the more fixed parameter values from input array "Phys_par"
dz = Phys_par(1); %grid stepsize (m)

zm = In_Z(end); %max depth
zz = [0:dz:zm-dz]'; %solution depth domain

Kz_K1 = Phys_par(2); % open water diffusion parameter (-)
Kz_K1_ice = Phys_par(3); % under ice diffusion parameter (-)
Kz_N0 = Phys_par(4); % min. stability frequency (s-2)
C_shelter = Phys_par(5); % wind shelter parameter (-)
lat = Phys_par(6); %latitude (decimal degrees)
lon = Phys_par(7); %longitude (decimal degrees)
alb_melt_ice = Phys_par(8);   %albedo of melting ice (-)
alb_melt_snow = Phys_par(9); %albedo of melting snow (-)
PAR_sat = Phys_par(10);         %PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
f_par = Phys_par(11);           %Fraction of PAR in incoming solar radiation (-)
beta_chl = Phys_par(12);        %Optical cross_section of chlorophyll (m2 mg-1)
lambda_i = Phys_par(13);       %PAR light attenuation coefficient for ice (m-1)
lambda_s = Phys_par(14);       %PAR light attenuation coefficient for snow (m-1)
F_sed_sld = Phys_par(15);      %volume fraction of solids in sediment (= 1-porosity)
I_scV = Phys_par(16); %scaling factor for inflow volume (-)
I_scT = Phys_par(17); %scaling coefficient for inflow temperature (-)
I_scC = Phys_par(18); %scaling factor for inflow concentration of C (-)
I_scS = Phys_par(19); %scaling factor for inflow concentration of S (-)
I_scTP = Phys_par(20); %scaling factor for inflow concentration of total P (-)
I_scDOP = Phys_par(21); %scaling factor for inflow concentration of diss. organic P (-)
I_scChl = Phys_par(22); %scaling factor for inflow concentration of Chl a (-)
I_scDOC = Phys_par(23); %scaling factor for inflow concentration of DOC  (-)


% Unpack the more site specific parameter values from input array "Bio_par"

swa_b0 = Bio_par(1); % non-PAR light atteneuation coeff. (m-1)
swa_b1 = Bio_par(2); % DOC-related specific PAR attenuation coeff. of water (m2 mg-1)
S_res_epi = Bio_par(3);      %Particle resuspension mass transfer coefficient, epilimnion (m day-1, dry)
S_res_hypo = Bio_par(4);     %Particle resuspension mass transfer coefficient, hypolimnion (m day-1, dry)
H_sed = Bio_par(5);          %height of active sediment layer (m, wet mass)
Psat_L = Bio_par(6);           %Half saturation parameter for Langmuir isotherm
Fmax_L = Bio_par(7);    %Scaling parameter for Langmuir isotherm !!!!!!!!!!!!

w_s = Bio_par(8);            %settling velocity for S (m day-1)
w_chl = Bio_par(9);          %settling velocity for Chl a (m day-1)
Y_cp = Bio_par(10);            %yield coefficient (chlorophyll to carbon) * (carbon to phosphorus) ratio (-)
m_twty = Bio_par(11);          %loss rate (1/day) at 20 deg C
g_twty = Bio_par(12);          %specific growth rate (1/day) at 20 deg C
k_twty = Bio_par(13);          %specific Chl a to P transformation rate (1/day) at 20 deg C
dop_twty = Bio_par(14);        %specific DOP to P transformation rate (day-1) at 20 deg C
P_half = Bio_par(15);          %Half saturation growth P level (mg/m3)

%NEW!!!===parameters for the 2 group of chlorophyll variable
PAR_sat_2 = Bio_par(16);        %PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
beta_chl_2 = Bio_par(17);       %Optical cross_section of chlorophyll (m2 mg-1)
w_chl_2 = Bio_par(18);          %Settling velocity for Chl a (m day-1)
m_twty_2 = Bio_par(19);         %Loss rate (1/day) at 20 deg C
g_twty_2 = Bio_par(20);         %Specific growth rate (1/day) at 20 deg C
P_half_2 = Bio_par(21);         %Half saturation growth P level (mg/m3)

%===========

% Parameters for C system
theta_DOC = Bio_par(22);        %Temperature adjustment coefficient for DOC degradation, T >= 10°C
theta_DOC_ice = Bio_par(23);    %Temperature adjustment coefficient for DOC degradation, T < 10°C
theta_sed = Bio_par(24);        %Temperature adjustment coefficient for sediment OC degradation, T >= 10°C
theta_sed_ice = Bio_par(25);    %Temperature adjustment coefficient for sediment OC degradation, T < 10°C
I_scDIC = Bio_par(26);          %Scaling factor for inflow concentration of DIC  (-)
I_scDO = Bio_par(27);           %Scaling factor for inflow concentration of dissolved oxygen  (-)
DOC_excr = Bio_par(28);         %DOC excretion fraction (-)
PP_resp = Bio_par(29);          %Phytoplankton respiration fraction (-)
C_per_P_alloc = Bio_par(30);    %C to P ratio in allochthonous particulate organic matter (mg C (mg P)-1)
C_P_org_sed = Bio_par(31);      %C to P ratio in dry nonliving organic sediment (mg C (mg P)-1)
k_POC_twty = Bio_par(32);       %POC degradation rate in sediment at 20°C (d-1)
k_POC1 = Bio_par(33);           %Labile (autochthonous) POC fragmentation rate at 20°C  (d-1)
k_POC2 = Bio_par(34);           %Semi-labile (allochthonous) POC fragmentation rate at 20°C  (d-1)
k_DOC1 = Bio_par(35);           %Labile (autochthonous) DOC degradation rate at 20°C (d-1)
k_DOC2 = Bio_par(36);           %Semi-labile DOC degradation rate at 20°C (d-1)
theta_cold = Bio_par(37);       %Temperature adjustment coefficient for OC degradation, T < 4°C (-)
PQuot = Bio_par(38);            %Photosynthetic quotient (mol O2 (mol CO2)-1)
RQuot = Bio_par(39);            %Respiratory quotient (mol CO2 (mol O2)-1)
CO2air = Bio_par(40);           %Atmospheric CO2 mixing ratio (ppm)
theta_POC_ice = theta_DOC_ice;  %Temperature adjustment coefficient for POC degradation, T < 10°C (-)
% ====== Other variables/parameters not read from the input file:

Nz=length(zz); %total number of layers in the water column
N_sed=26; %total number of layers in the sediment column

theta_m = exp(0.1*log(2));    %Chl loss and growth rate & POC fragmentation (T >= 10°C) rate parameter base, ~1.072
e_par = 240800;               %Average energy of PAR photons (J mol-1)

if(DOC_attenuation_switch == 0)
    swa_b1 = swa_b1*10*1000;
end

% diffusion parameterisation exponents
Kz_b1 = 0.43;
Kz_b1_ice  = 0.43;

% ice & snow parameter values
rho_fw=1000;        %density of freshwater (kg m-3)
rho_ice=910;        %ice (incl. snow ice) density (kg m-3)
rho_new_snow=250;   %new-snow density (kg m-3)
max_rho_snow=550;   %maximum snow density (kg m-3)
L_ice=333500;       %latent heat of freezing (J kg-1)
K_ice=2.1;          %ice heat conduction coefficient (W m-1 K-1)
C1=7.0;             %snow compaction coefficient #1
C2=21.0;            %snow compaction coefficient #2

Tf=0;               %water freezing point temperature (deg C)

F_OM=1e+6*0.010;    %mass fraction [mg kg-1] of P of dry organic matter (assuming 50% of C, and Redfield ratio)

K_sed=0.015;      %thermal diffusivity of the sediments (m2 day-1) 0.035
rho_sed=2500;      %bulk density of the inorganic solids in sediments (kg m-3)
rho_org=1000;      %bulk density of the organic solids in sediments (kg m-3)
cp_sed=1500;       %specific heat capacity of the sediments (J kg-1 K-1) 1000

ksw=1e-3; %sediment pore water mass transfer coefficient (m/d)
Fmax_L_sed=Fmax_L;
Fstable=655; % Inactive P conc. in inorg. particles (mg/kg dw);

Frazil2Ice_tresh=0.013;  % treshold (m) where frazil is assumed to turn into a solid ice cover NEW!!!

flocc = 0.0019; %DOC flocculation rate
F_OC = 0.5e6; %C conc. in nonliving organic sediment particles (mg kg-1 dry w.)
C_per_Chla = 50/Y_cp; % C to Chl a ratio in phytoplankton (mg C (mg Chl a)-1)
C_per_P_pp = 50; % C to P ratio in autochthonous particulate organic matter (mg C (mg P)-1)
%=======

% Allocate and initialise output data matrices
Qst = zeros(3,length(tt));
Kzt = zeros(Nz,length(tt));
Tzt = zeros(Nz,length(tt));
Czt = zeros(Nz,length(tt));
Szt = zeros(Nz,length(tt));
Pzt = zeros(Nz,length(tt));
Chlzt = zeros(Nz,length(tt));
PPzt = zeros(Nz,length(tt));
DOPzt = zeros(Nz,length(tt));
DOCzt = zeros(Nz,length(tt));
DICzt = zeros(Nz,length(tt));
CO2zt = zeros(Nz,length(tt));
O2zt = zeros(Nz,length(tt));
POCzt = zeros(Nz,length(tt));
Salzt = zeros(Nz,length(tt));
O2_sat_relt = zeros(Nz,length(tt));
O2_sat_abst = zeros(Nz,length(tt));
T_sedt = zeros(Nz,length(tt),4);
Qzt_sed = zeros(Nz,length(tt));
lambdazt = zeros(Nz,length(tt));
P3zt_sed = zeros(Nz,length(tt),11); %3-D
P3zt_sed_sc = zeros(Nz,length(tt),3); %3-D
His = zeros(8,length(tt)); %NEW!!!
MixStat = zeros(23,length(tt));
% Fokema
DOCzt1=zeros(Nz,length(tt)); %Fokema-model subpool 1
DOCzt2=zeros(Nz,length(tt)); %Fokema-model subpool 2
DOCzt3=zeros(Nz,length(tt)); %Fokema-model subpool 3
DOCtfrac=zeros(Nz,length(tt),3); %Fokema-model subpool fractions
Daily_BB1t=zeros(Nz,length(tt)); %Fokema-model subpool 1 daily bacterial decomposition
Daily_BB2t=zeros(Nz,length(tt)); %Fokema-model subpool 2 daily bacterial decomposition
Daily_BB3t=zeros(Nz,length(tt)); %Fokema-model subpool 3 daily bacterial decomposition
Daily_PBt=zeros(Nz,length(tt)); %Fokema-model daily photobleaching

surfaceflux = zeros(1,length(tt)); %CO2 surface flux
CO2_eqt = zeros(1,length(tt));     %CO2 equilibrium concentration
CO2_ppmt = zeros(1,length(tt));    %CO2 fraction in air
K0t = zeros(1,length(tt));         %CO2 solubility coefficient

O2fluxt = zeros(1,length(tt));     %oxygen surface flux
O2_eqt = zeros(1,length(tt));      %O2 equilibrium concentration
K0_O2t = zeros(1,length(tt));      %O2 solubility coefficient
dO2Chlt = zeros(Nz,length(tt));    %Oxygen change due to phytoplankton (mg m-3))
POC1tfrac = zeros(Nz,length(tt));    %Oxygen consumption due to BOD (mg m-3))
dO2SODt = zeros(Nz,length(tt));    %Oxygen consumption due to SOD (mg m-3))
dO2DOCt = zeros(Nz,length(tt));    %Oxygen consumption due to DOC (mg m-3))
pHt = zeros(Nz,length(tt));        %pH
testi3t = zeros(Nz,length(tt),6); %3-D
% Initial profiles

Az = interp1(In_Z,In_Az,zz);
Vz = dz * (Az + [Az(2:end); 0]) / 2;

T0 = interp1(In_Z,In_Tz,zz+dz/2); % Initial temperature distribution (deg C)
C0 = interp1(In_Z,In_Cz,zz+dz/2); % Initial  chlorophyll (group 2) distribution (mg m-3)
S0 = interp1(In_Z,In_Sz,zz+dz/2); % Initial passive sedimenting tracer (or suspended inorganic matter) distribution (kg m-3)
TP0 = interp1(In_Z,In_TPz,zz+dz/2);	% Initial total P distribution (incl. DOP & Chla & Cz) (mg m-3)
DOP0 = interp1(In_Z,In_DOPz,zz+dz/2);	% Initial dissolved organic P distribution (mg m-3)
Chl0 = interp1(In_Z,In_Chlz,zz+dz/2);	% Initial chlorophyll (group 2) distribution (mg m-3)
DOC0 = interp1(In_Z,In_DOCz,zz+dz/2);	% Initial DOC distribution (mg m-3)
DIC0 = interp1(In_Z,In_DICz,zz+dz/2);   % Initial DIC distribution (mg m-3)
O20 = interp1(In_Z,In_O2z,zz+dz/2);   % Initial oxygen distribution (mg m-3)
POC0 = interp1(In_Z,In_POCz,zz+dz/2);	% Initial POC distribution (mg m-3)
TP0_sed = interp1(In_Z,In_TPz_sed,zz+dz/2); % Initial total P distribution in bulk wet sediment ((mg m-3); particles + porewater)
FNLOM0 = interp1(In_Z,In_FNLOM,zz+dz/2);     % Initial sediment solids volume fraction of nonliving organic matter (-)
FIM0 = interp1(In_Z,In_FIM,zz+dz/2);     % Initial sediment solids volume fraction of inorganic matter (-)
pH0 = interp1(In_Z,In_pH,zz+dz/2);     % Initial sediment solids volume fraction of inorganic matter (-)
Hplus0 = 10.^(-pH0); % Initial pH (-)
Sal0 = interp1(In_Z,0*In_Tz,zz+dz/2);

VolFrac=1./(1+(1-F_sed_sld)./(F_sed_sld*FIM0)); %volume fraction: inorg sed. / (inorg.sed + pore water)
%keyboard
if (any(FIM0<0)||any(FIM0>1))
    error('Initial fraction of inorganic matter in sediments must be between 0 and 1')
end

if (any(ksw>(H_sed*(1-F_sed_sld))))
    error('Parameter ksw is larger than the volume (thickness) of porewater')
end  %OBS! Ideally should also be that the daily diffused porewater should not be larger
%than the corresponding water layer volume, but this seems very unlike in practise

Tz = T0;
Cz = C0; % (mg m-3)
Sz = S0; % (kg m-3)
Chlz = Chl0;  % (mg m-3)
DOPz = DOP0;  % (mg m-3)
DOCz = DOC0;   % (mg m-3)
DICz = DIC0;   % (mg m-3)
O2z = O20;   % (mg m-3)
POCz = POC0;   % (mg m-3)
Salz = Sal0;
F_IM = FIM0; %initial VOLUME fraction of inorganic particles of total dry sediment solids
F_NLOM = FNLOM0; % initial VOLUME fraction of nonliving organic particles of total dry sediment solids
F_LOM = 1 - (F_IM + F_NLOM); % initial VOLUME fraction of living organic particles of total dry sediment solids
Hplusz = Hplus0;

theta_sod_factor_thresh = theta_sed_ice^(4-20);
theta_POC_factor_thresh = theta_POC_ice^(4-20);

POCz1 = OC_frac0(3)*POCz; %Initial autochthonous POC fraction in water column
POCz2 = POCz-POCz1; %Initial allochthonous POC fraction in water column

DOC1_frac = OC_frac0(1); %Initial autochthonous DOC fraction in water column
DOC2_frac = OC_frac0(2); %Initial semi-labile allochthonous DOC fraction in water column
DOC3_frac = 1-(DOC1_frac+DOC2_frac); %Initial semi-labile refractory DOC fraction in water column

DOCz1 = DOC1_frac.*DOCz; %DOC subpools
DOCz2 = DOC2_frac.*DOCz;
DOCz3 = DOC3_frac.*DOCz;

POP0 = POCz1./(C_per_P_pp) + POCz2./C_per_P_alloc;

[Pz, trash] =Ppart(S0./rho_sed,TP0-DOP0-POP0-((Chl0 + C0)./Y_cp),Psat_L,Fmax_L,rho_sed,Fstable);  % (mg m-3) NEW!!!
PPz = TP0-DOP0-POP0-((Chl0 + C0)./Y_cp)-Pz; % (mg m-3) NEW!!!
%keyboard
if any((TP0-DOP0-POP0-(Chl0 + C0)./Y_cp-S0*Fstable)<0) %NEW!!!
    keyboard
    error('Sum of initial DOP, stably particle bound P, and P contained in Chl (both groups) a cannot be larger than TP')
end

Chlsz_store = 12500*ones(Nz,1); %Chlsz_store: Chla conc. in organic sediment particles (mg kg-1 dry w.)
%Chlsz_store = Chl0_sed./(rho_org*F_sed_sld*F_LOM); %(mg kg-1 dry w.)

C_P_org_sed = C_P_org_sed*ones(Nz,1);

POP0_sed = rho_org*F_sed_sld.*F_NLOM*F_OC./C_P_org_sed; %Sediment nonliving particulate organic phosphorus in bulk wet sediment
Chl0_sed = Chlsz_store*rho_org*F_sed_sld.*F_LOM; % Initial chlorophyll (group 1+2) distribution in bulk wet sediment (mg m-3)

%Tähän sedimentin päivittäinen POP kuivapitoisuutena
POPz_sed = POP0_sed./F_sed_sld; %(mg P) (m-3 dry sed.)

if any((TP0_sed-DOP0-POP0_sed-Chl0_sed./Y_cp-VolFrac*rho_sed*Fstable)<0)
    keyboard
end

if any((TP0_sed-DOP0-POP0_sed-Chl0_sed./Y_cp-VolFrac*rho_sed*Fstable)<0)
    error('Sum of initial DOP stably, particle bound P, and P contained in Chl_sed a cannot be larger than TP_sed')
end

%== P-partitioning in sediments==
%Pdz_store: %diss. inorg. P in sediment pore water (mg m-3)
%Psz_store: %P conc. in inorganic sediment particles (mg kg-1 dry w.)

[Pdz_store, Psz_store]=Ppart(VolFrac,TP0_sed-POP0_sed-(Chl0_sed./Y_cp)-DOP0,Psat_L,Fmax_L_sed,rho_sed,Fstable);

%T_sed_mdl = [10 10 10 10 10 10 9 9 8 8 8 7 7 7 6 6 6 5.5 5.5 5 4.5 4.5 4.5 4.5 4 4 4 4];
%T_sed_btm = [8.1 9.1 9.4 9.4 9.2 8.7 8.1 7.6 7.4 7.1 6.4 5.9 5.6 5.2 5.0 4.8 4.7 4.6 4.5 4.4 4.3 4.2 4.1 4.1 4.1 4.1 4.1 4.1];
T_sed_btm = [8.3 8.6 8.6 8.6 8.6 8.3 8.2 8.0 7.9 7.8 7.5 7.2 7.0 6.7 6.5 6.2 6.0 5.9 5.8 5.8 5.6 5.6 5.6 5.7 5.7 5.7 5.7 5.7];
T_sed_mdl = 0.146*(T_sed_btm-Tz')+T_sed_btm;
% assume linear initial temperature profile in sediment (4 deg C at the bottom)
clear Tzy_sed
for j=1:Nz
    Tzy_sed(:,j) = interp1([0.2 2 3.5 6 10], [Tz(j) T_sed_btm(j) T_sed_mdl(j) T_sed_btm(j) T_sed_btm(j)], [0.2:0.2:2 2.5:0.5:10])';
    %Tzy_sed(:,j) = interp1([0.2 2 10], [Tz(j) 6 8], [0.2:0.2:2 2.5:0.5:10])';
end

S_resusp=S_res_hypo*ones(Nz,1); %hypolimnion resuspension assumed on the first time step

rho_snow=rho_new_snow;   %initial snow density (kg m-3)
Tice=NaN;                %ice surface temperature (initial value, deg C)
XE_melt=0;               %energy flux that is left from last ice melting (initial value, W m-2)
XE_surf=0;               %energy flux from water to ice (initial value,  J m-2 per day)

%Initialisation of ice & snow variables
Hi=Ice0(1);               %total ice thickness (initial value, m)
WEQs=(rho_snow/rho_fw)*Ice0(2); %snow water equivalent  (initial value, m)
Hsi=0;                %snow ice thickness (initial value = 0 m)
HFrazil=0;              % (initial value, m) NEW!!!


if ((Hi<=0)&&(WEQs>0))
    error('Mismatch in initial ice and snow thicknesses')
end

if (Hi<=0)
    IceIndicator=0;     %IceIndicator==0 means no ice cover
else
    IceIndicator=1;
end

pp=1; %initial indexes for ice freezing/melting date arrays
qq=1;
DoF=[]; %initialize
DoM=[]; %initialize

% >>>>>> Start of the time loop >>>>>>
Resuspension_counter=zeros(Nz,1); %kg
Sedimentation_counter=zeros(Nz,1); %kg
POC_sedim_counter=zeros(Nz,1); %kg
POC_resusp_counter=zeros(Nz,1); %kg
SS_decr=0; %kg
alb_melt_snow_0 = alb_melt_snow;
naytto = -122;
DICpaiva = -112;

for i = 1:length(tt)
    if(datenum(M_start)== datenum([2013 1 8]) && i==367);Hi=0.01;Hsi=0;WEQs=0;end
    if(rho_snow>=0.9*max_rho_snow)
        alb_melt_snow = alb_melt_snow_0/2;
    else
        alb_melt_snow = alb_melt_snow_0;
    end
    
    % Surface heat fluxes (W m-2), wind stress (N m-2) & daylight fraction (-), based on Air-Sea Toolbox
    [Qsw,Qlw,Qsl,tau,DayFrac,DayFracHeating,Qs,Ql] = heatflux_v12(tt(i),Wt(i,1),Wt(i,2),Wt(i,3),Wt(i,4),Wt(i,5),Wt(i,6),Tz(1), ...
        lat,lon,WEQs,Hi,alb_melt_ice,alb_melt_snow,albedot1);     %Qlw and Qsl are functions of Tz(1)
    
    % Calculate total mean PAR and non-PAR light extinction coefficient in water (background + due to Chl a)
    lambdaz_wtot_avg=zeros(Nz,1);
    lambdaz_NP_wtot_avg=zeros(Nz,1);
    
    if(DOC_attenuation_switch == 1)
        DOCz = DOCz1 + DOCz2 + DOCz3;
        if (selfshading_switch==1)
            lambdaz_wtot=swa_b1*DOCz + beta_chl*Chlz + beta_chl_2*Cz; %at layer z
            %keyboard
            %lambdaz_NP_wtot=swa_b0 * ones(Nz,1) + beta_chl*Chlz + beta_chl_2*Cz; %at layer z
            for j=1:Nz
                lambdaz_wtot_avg(j)=mean(swa_b1*DOCz(1:j) + beta_chl*Chlz(1:j) + beta_chl_2*Cz(1:j)); %average down to layer z
                lambdaz_NP_wtot_avg(j)=mean(swa_b0 * ones(j,1) + beta_chl*Chlz(1:j) + beta_chl_2*Cz(1:j)); %average down to layer z
            end
        else %no Chl selfshading effect
            lambdaz_wtot=swa_b1*DOCz;
            %lambdaz_NP_wtot=swa_b0 * ones(Nz,1);
            for j=1:Nz
                lambdaz_wtot_avg(j)=mean(swa_b1*DOCz(1:j)); %average down to layer z
            end
            lambdaz_NP_wtot_avg=swa_b0 * ones(Nz,1);
        end %if selfshading...
        
    else
        
        %NEW!!! below additional term for chlorophyll group 2
        if (selfshading_switch==1)
            lambdaz_wtot=swa_b1 * ones(Nz,1) + beta_chl*Chlz + beta_chl_2*Cz; %at layer z
            lambdaz_NP_wtot=swa_b0 * ones(Nz,1) + beta_chl*Chlz + beta_chl_2*Cz; %at layer z
            for j=1:Nz
                lambdaz_wtot_avg(j)=mean(swa_b1 * ones(j,1) + beta_chl*Chlz(1:j) + beta_chl_2*Cz(1:j)); %average down to layer z
                lambdaz_NP_wtot_avg(j)=mean(swa_b0 * ones(j,1) + beta_chl*Chlz(1:j) + beta_chl_2*Cz(1:j)); %average down to layer z
            end
        else %constant with depth
            lambdaz_wtot=swa_b1 * ones(Nz,1);
            lambdaz_wtot_avg=swa_b1 * ones(Nz,1);
            lambdaz_NP_wtot=swa_b0 * ones(Nz,1);
            lambdaz_NP_wtot_avg=swa_b0 * ones(Nz,1);
        end %if selfshading...
        
    end %if DOC_attenuation...
    
    if(IceIndicator==0)
        IceSnowAttCoeff=1; %no extra light attenuation due to snow and ice
    else    %extra light attenuation due to ice and snow
        IceSnowAttCoeff=exp(-lambda_i * (Hi)) * exp(-lambda_s * ((rho_fw/rho_snow)*WEQs));
        SnowAttCoeff=exp(-lambda_s * (rho_fw/rho_snow)*WEQs);
        %IceAttCoeff=exp(-lambda_i * 0.5*Hi) * exp(-lambda_s * (rho_fw/rho_snow)*WEQs);
    end
    T_testi = Tz;
    Tprof_prev=Tz; %temperature profile at previous time step (for convection_v12_1a.m)
    if(i==naytto);display(['Tprof_prev = ',num2str(Tprof_prev')]);end
    if(density_switch == 0)
        rho = polyval(ies80,max(0,Tz(:))) + min(Tz(:),0);  % Density (kg/m3)
    else
        rho = sw_dens(Salz(:),max(0,Tz(:)),0) + min(Tz(:),0);
    end
    
    % Sediment vertical heat flux, Q_sed
    % (averaged over the whole top area of the layer, although actually coming only from the "sides")
    if (sediment_heatflux_switch==1)
        % update top sediment temperatures
        dz_sf = 0.2; %fixed distance between the two topmost sediment layers (m)
        Tzy_sed(1,:) = Tz';
        Tzy_sed_upd = sedimentheat_v11(Tzy_sed, K_sed, dt);
        Tzy_sed=Tzy_sed_upd;
        Qz_sed=K_sed*rho_sed*cp_sed*(1/dz_sf)*(-diff([Az; 0])./Az) .* (Tzy_sed(2,:)'-Tzy_sed(1,:)'); %(J day-1 m-2)
        T_sedim = Tzy_sed(10,:)';
        T_sedim_m = Tzy_sed(18,:)';
        T_sedim_b = Tzy_sed(end,:)';
        T_sedim_t = Tzy_sed(4,:)';
        %positive heat flux => from sediment to water
    else
        Qz_sed = zeros(Nz,1);
    end
    
    Cw = 4.18e+6;	% Volumetric heat capacity of water (J K-1 m-3)
    
    %Heat sources/sinks:
    %Total attenuation coefficient profile, two-band extinction, PAR & non-PAR
    Par_Attn=exp([0; -lambdaz_wtot_avg] .* [zz; zz(end)+dz]);
    NonPar_Attn=exp([0; -lambdaz_NP_wtot_avg] .* [zz; zz(end)+dz]);
    
    Attn_z=(-f_par * diff([1; ([Az(2:end);0]./Az).*Par_Attn(2:end)]) + ...
        (-(1-f_par)) * diff([1; ([Az(2:end);0]./Az).*NonPar_Attn(2:end)])); %NEW (corrected 210807)
    
    if(IceIndicator==0)
        % 1) Vertical heating profile for open water periods (during daytime heating)
        Qz = (Qsw + XE_melt) * Attn_z; %(W m-2)
        Qz(1) = Qz(1) + DayFracHeating*(Qlw + Qsl); %surface layer heating
        XE_melt=0; %Reset
        dT = Az .* ((60*60*24*dt) * Qz + DayFracHeating*Qz_sed) ./ (Cw * Vz); %Heat source (K day-1) (daytime heating, ice melt, sediment);
        
        % === Frazil ice melting, NEW!!! === %
        postemp=find(dT>0);
        if (isempty(postemp)==0)
            RelT=dT(postemp)./sum(dT(postemp));
            HFrazilnew=max(0, HFrazil - sum(dT(postemp))*1/((Az(1)*rho_ice*L_ice)/(Cw * Vz(1)))); %
            sumdTnew = max(0, sum(dT(postemp))-(HFrazil*Az(1)*rho_ice*L_ice)/(Cw * Vz(1)));
            dT(postemp)=RelT.*sumdTnew;
            HFrazil=HFrazilnew;
        end
        % === === ===
    else
        % Vertical heating profile for ice-covered periods (both day- and nighttime)
        Qz = Qsw * IceSnowAttCoeff * Attn_z; %(W/m2)
        dT = Az .* ((60*60*24*dt) * Qz + Qz_sed) ./ (Cw * Vz); %Heat source (K day-1) (solar rad., sediment);
    end
    
    [CO2z,~] = carbonequilibrium(DICz,Tz,-log10(Hplusz));
    Tz_IC_old = Tz; %old temperature for IC system equilibrium calculation
    if(i==DICpaiva);keyboard;end
    Tz = Tz + dT;        %Temperature change after daytime surface heatfluxes (or whole day in ice covered period)
    if(i==naytto);display(['Päivävuo ja sedimentti ',num2str(Tz')]);end
    
    [CO2z,DICz,Hplusz] = DICsystem_new(DICz,CO2z,Tz_IC_old,Tz,Hplusz,i);
    if(i==DICpaiva);keyboard;end
    % Convective mixing adjustment (mix successive layers until stable density profile)
    % and
    % Spring/autumn turnover (don't allow temperature jumps over temperature of maximum density)
    [Tz,Cz,Sz,Pz,Chlz,PPz,DOPz,DOCz1,DOCz2,DOCz3,DICz,O2z,POCz1,POCz2,Hplusz,CO2z] = convection_v12_1aIC(Tz,Cz,Sz,Pz,Chlz,PPz,DOPz,DOCz1,DOCz2,DOCz3,DICz,O2z,POCz1,POCz2,Hplusz,CO2z,Tprof_prev,Vz,Cw,f_par,lambdaz_wtot_avg,zz,swa_b0,tracer_switch,1,i);
    
    if(i==DICpaiva);keyboard;end
    Tprof_prev=Tz; %NEW!!! Update Tprof_prev
    if(i==naytto);display(['Konvektio ',num2str(Tz')]);end
    if(IceIndicator==0)
        % 2) Vertical heating profile for open water periods (during nighttime heating)
        [Qsw,Qlw_2,Qsl_2,tau,DayFrac,DayFracHeating,Qs_2,Ql_2] = heatflux_v12(tt(i),Wt(i,1),Wt(i,2),Wt(i,3),Wt(i,4),Wt(i,5),Wt(i,6),0.5*(2*T_testi(1)+0*Tz(1)), ...
            lat,lon,WEQs,Hi,alb_melt_ice,alb_melt_snow,albedot1); %Qlw and Qsl are functions of Tz(1)
        Qz(1) = (1-DayFracHeating)*(Qlw_2 + Qsl_2); %surface layer heating
        Qz(2:end)=0; %No other heating below surface layer
        dT = Az .* ((60*60*24*dt) * Qz + (1-DayFracHeating)*Qz_sed) ./ (Cw * Vz); %Heat source (K day-1) (longwave & turbulent fluxes);
        %keyboard
        % === NEW!!! frazil ice melting === %
        postemp=find(dT>0);
        if (isempty(postemp)==0)
            %disp(['NOTE: positive night heat flux at T=' num2str(Tz(postemp),2)]) %NEW
            RelT=dT(postemp)./sum(dT(postemp));
            HFrazilnew=max(0, HFrazil - sum(dT(postemp))*1/((Az(1)*rho_ice*L_ice)/(Cw * Vz(1)))); %
            sumdTnew = max(0, sum(dT(postemp))-(HFrazil*Az(1)*rho_ice*L_ice)/(Cw * Vz(1)));
            dT(postemp)=RelT.*sumdTnew;
            HFrazil=HFrazilnew;
        end
        % === === ===
        
        Tz_IC_old = Tz;
        Tz = Tz + dT;         %Temperature change after nighttime surface heatfluxes
        [CO2z,DICz,Hplusz] = DICsystem_new(DICz,CO2z,Tz_IC_old,Tz,Hplusz,i);
        if(i==naytto);display(['Yövuo ',num2str(Tz')]);end
        % Convective mixing adjustment (mix successive layers until stable density profile)
        % and
        % Spring/autumn turnover (don't allow temperature jumps over temperature of maximum density)
        [Tz,Cz,Sz,Pz,Chlz,PPz,DOPz,DOCz1,DOCz2,DOCz3,DICz,O2z,POCz1,POCz2,Hplusz,CO2z] = convection_v12_1aIC(Tz,Cz,Sz,Pz,Chlz,PPz,DOPz,DOCz1,DOCz2,DOCz3,DICz,O2z,POCz1,POCz2,Hplusz,CO2z,Tprof_prev,Vz,Cw,f_par,lambdaz_wtot_avg,zz,swa_b0,tracer_switch,1,i);
        if(i==naytto);display(['Konvektio ',num2str(Tz')]);end
        if(i==DICpaiva);keyboard;end
        Qlw = DayFracHeating*Qlw + (1-DayFracHeating)*Qlw_2; %total amounts, only for output purposes
        Qsl = DayFracHeating*Qsl + (1-DayFracHeating)*Qsl_2; %total amounts, only for output purposes
        Qs = DayFracHeating*Qs + (1-DayFracHeating)*Qs_2;
        Ql = DayFracHeating*Ql + (1-DayFracHeating)*Ql_2;
    end
    %keyboard
    % Vertical turbulent diffusion
    g   = 9.81;							% Gravity acceleration (m s-2)
    if(density_switch == 0)
        rho = polyval(ies80,max(0,Tz(:))) + min(Tz(:),0);  % Water density (kg m-3)
    else
        rho = sw_dens(Salz,max(0,Tz(:)),0) + min(Tz(:),0);
    end
    % Note: in equations of rho it is assumed that every supercooled degree lowers density by
    % 1 kg m-3 due to frazil ice formation (probably no practical meaning, but included for "safety")
    
    N2  = g * (diff(log(rho)) ./ diff(zz));	% Brunt-Vaisala frequency (s-2) for level (zz+1)
    if (IceIndicator==0)
        Kz  = Kz_K1 * max(Kz_N0, N2).^(-Kz_b1);	% Vertical diffusion coeff. in ice free season (m2 day-1)
        Kz_O2 = Kz; %Vertical diffusion coeff. in ice free season for oxygen (m2 day-1) == NEW NEW NEW==
        % for level (zz+1)
    else
        Kz  = Kz_K1_ice * max(Kz_N0, N2).^(-Kz_b1_ice); % Vertical diffusion coeff. under ice cover (m2 day-1)
        Kz_O2 = Kz; %Vertical diffusion coeff. under ice cover for oxygen (m2 day-1) == NEW NEW NEW==
        % for level (zz+1)
    end
    
    Tz_IC_old = Tz;
    
    Fi = tridiag_DIF_v11([NaN; Kz],Vz,Az,dz,dt); %Tridiagonal matrix for general diffusion
    
    Tz = Fi \ (Tz);        %Solving new temperature profile (diffusion, sources/sinks already added to Tz above)
    if(i==naytto);display(['Diffuusio ',num2str(Tz')]);end
    
    [CO2z,DICz,Hplusz] = DICsystem_new(DICz,CO2z,Tz_IC_old,Tz,Hplusz,i);
    
    if(i==DICpaiva);keyboard;end
    % Convective mixing adjustment (mix successive layers until stable density profile)
    % (don't allow temperature jumps over temperature of maximum density, no summer/autumn turnover here!)
    [Tz,Cz,Sz,Pz,Chlz,PPz,DOPz,DOCz1,DOCz2,DOCz3,DICz,O2z,POCz1,POCz2,Hplusz,CO2z] = convection_v12_1aIC(Tz,Cz,Sz,Pz,Chlz,PPz,DOPz,DOCz1,DOCz2,DOCz3,DICz,O2z,POCz1,POCz2,Hplusz,CO2z,Tprof_prev,Vz,Cw,f_par,lambdaz_wtot_avg,zz,swa_b0,tracer_switch,0,i);
    if(i==naytto);display(['Konvektio ',num2str(Tz')]);end
    
    if(i==DICpaiva);keyboard;end
    dSz_porg = rho_org*S_resusp.*F_NLOM.*(-diff([Az; 0])./Vz);  %Dry organic particle resuspension source from sediment (kg m-3 day-1)
    dPOC_res = dSz_porg.*F_OC;  %POC resuspension source from sediment resusp. ((kg m-3 day-1)*(mg kg-1) = mg m-3 day-1);
    
    %POC degradation to DOC
    POC_t_fctr = NaN*ones(length(zz),1);
    for k = 1:length(Tz)
        if(Tz(k)<10)
            if(Tz(k)<4)
                POC_t_fctr(k) = theta_POC_factor_thresh.*theta_cold.^(max(0,Tz(k))-4);
            else
                POC_t_fctr(k) = theta_POC_ice^(Tz(k)-20);
            end
        else
            POC_t_fctr(k) = theta_m^(Tz(k)-20);
        end
    end
    
    POC_sed_t_fctr = NaN*ones(length(zz),1);
    for k = 1:length(Tz)
        if(Tz(k)<10)
            if(Tz(k)<4)
                POC_sed_t_fctr(k) = theta_sod_factor_thresh.*theta_cold.^(max(0,Tz(k))-4);
            else%tähän ero veden POC-hajoamiseen
                POC_sed_t_fctr(k) = theta_sed_ice^(Tz(k)-20);
            end
        else
            POC_sed_t_fctr(k) = theta_sed^(Tz(k)-20);
        end
    end
    
    dPOC1 = k_POC1*POCz1.*POC_t_fctr;
    dPOC2 = k_POC2*POCz2.*POC_t_fctr;
    POCz1 = POCz1 - dPOC1;
    POCz2 = POCz2 - dPOC2;
    
    dDOP =  dop_twty * DOPz .* theta_m.^(Tz-20);  %Mineralisation to P
    DOPz = Fi \ (DOPz - dDOP);         %Solving new dissolved inorganic P profile (diffusion)
    
    %Dissolved inorganic phosphorus
    dP = dDOP + dPOC1./(C_per_P_pp) + dPOC2./C_per_P_alloc;
    Pz = Pz + dP; %Solving new dissolved inorganic P profile (diffusion)
    
    % NEW!!! === Code rearranging
    % Calculate again the total mean PAR light extinction coefficient in water (background + due to Chl a)
    lambdaz_wtot_avg=zeros(Nz,1);
    
    if(DOC_attenuation_switch == 1)
        
        %NEW!!! below additional term for chlorophyll group 2
        if (selfshading_switch==1)
            lambdaz_wtot=swa_b1*DOCz + beta_chl*Chlz + beta_chl_2*Cz; %at layer z.
            for j=1:Nz
                lambdaz_wtot_avg(j)=mean(swa_b1*DOCz(1:j) + beta_chl*Chlz(1:j) + beta_chl_2*Cz(1:j)); %average down to layer z
            end
        else %constant with depth
            lambdaz_wtot=swa_b1*DOCz;
            for j=1:Nz
                lambdaz_wtot_avg(j)=mean(swa_b1*DOCz(1:j)); %average down to layer z
            end
        end %if selfshading...
        
    else
        
        %NEW!!! below additional term for chlorophyll group 2
        if (selfshading_switch==1)
            lambdaz_wtot=swa_b1 * ones(Nz,1) + beta_chl*Chlz + beta_chl_2*Cz; %at layer z.
            for j=1:Nz
                lambdaz_wtot_avg(j)=mean(swa_b1 * ones(j,1) + beta_chl*Chlz(1:j) + beta_chl_2*Cz(1:j)); %average down to layer z
            end
        else %constant with depth
            lambdaz_wtot=swa_b1 * ones(Nz,1);
            lambdaz_wtot_avg=swa_b1 * ones(Nz,1);
        end %if selfshading...
        
    end
    
    %Photosynthetically Active Radiation (for chlorophyll group 1)
    H_sw_z=NaN*zeros(Nz,1);
    
    % ===== NEW!!! bug (when Dayfrac==0) fixed 071107
    if ((IceIndicator==0)&&(DayFrac>0))
        PAR_z=((3/2) / (e_par * DayFrac)) * f_par * Qsw  * exp(-lambdaz_wtot_avg .* zz);
        %Irradiance at noon (mol m-2 s-1) at levels zz
    elseif ((IceIndicator==1)&&(DayFrac>0))    %extra light attenuation due to ice and snow
        PAR_z=((3/2) / (e_par * DayFrac)) * IceSnowAttCoeff * f_par *...
            Qsw  * exp(-lambdaz_wtot_avg .* zz);
    else PAR_z=zeros(Nz,1); %DayFrac==0, polar night
    end
    % =====
    
    U_sw_z=PAR_z./PAR_sat; %scaled irradiance at levels zz
    inx_u=find(U_sw_z<=1); %undersaturated
    inx_s=find(U_sw_z>1);  %saturated
    
    H_sw_z(inx_u)=(2/3)*U_sw_z(inx_u);  %undersaturated
    
    dum_a=sqrt(U_sw_z);
    dum_b=sqrt(U_sw_z-1);
    H_sw_z(inx_s)=(2/3)*U_sw_z(inx_s) + log((dum_a(inx_s) + dum_b(inx_s))./(dum_a(inx_s) ...  %saturated
        - dum_b(inx_s))) - (2/3)*(U_sw_z(inx_s)+2).*(dum_b(inx_s)./dum_a(inx_s));
    
    
    %NEW!!!! modified for chlorophyll group 1
    Growth_bioz=g_twty*theta_m.^(Tz-20) .* (Pz./(P_half+Pz)) .* (DayFrac./(dz*lambdaz_wtot)) .* diff([-H_sw_z; 0]);
    Loss_bioz=m_twty*theta_m.^(Tz-20);
    R_bioz = Growth_bioz-Loss_bioz;
    
    %Photosynthetically Active Radiation (for chlorophyll group 2) NEW!!!
    H_sw_z=NaN*zeros(Nz,1);
    
    U_sw_z=PAR_z./PAR_sat_2; %scaled irradiance at levels zz
    inx_u=find(U_sw_z<=1); %undersaturated
    inx_s=find(U_sw_z>1);  %saturated
    
    H_sw_z(inx_u)=(2/3)*U_sw_z(inx_u);  %undersaturated
    
    dum_a=sqrt(U_sw_z);
    dum_b=sqrt(U_sw_z-1);
    H_sw_z(inx_s)=(2/3)*U_sw_z(inx_s) + log((dum_a(inx_s) + dum_b(inx_s))./(dum_a(inx_s) ...  %saturated
        - dum_b(inx_s))) - (2/3)*(U_sw_z(inx_s)+2).*(dum_b(inx_s)./dum_a(inx_s));
    
    Growth_bioz_2=g_twty_2*theta_m.^(Tz-20) .* (Pz./(P_half_2+Pz)) .* (DayFrac./(dz*lambdaz_wtot)) .* diff([-H_sw_z; 0]);
    Loss_bioz_2=m_twty_2*theta_m.^(Tz-20);
    R_bioz_2 = Growth_bioz_2-Loss_bioz_2;
    exinx = find( (Growth_bioz.*Chlz*dt + Growth_bioz_2.*Cz*dt)>(Y_cp*Pz) );
    
    if (isempty(exinx)==0)
        Growth_bioz_ratio = (Growth_bioz(exinx).*Chlz(exinx)*dt)./((Growth_bioz(exinx).*Chlz(exinx)*dt) + (Growth_bioz_2(exinx).*Cz(exinx)*dt)); %fraction of Growth rate 1 of total growth rate
        Growth_bioz(exinx) = Growth_bioz_ratio.*(Y_cp*Pz(exinx)./(Chlz(exinx)*dt));
        Growth_bioz_2(exinx) = (1-Growth_bioz_ratio).*(Y_cp*Pz(exinx)./(Cz(exinx)*dt));
        R_bioz = Growth_bioz-Loss_bioz;
        R_bioz_2 = Growth_bioz_2-Loss_bioz_2;
    end
    %================================
    
    if(i==DICpaiva);keyboard;end
    % Suspended solids, particulate inorganic P
    Fi_ad = tridiag_HAD_v11([NaN; Kz],w_s,Vz,Az,dz,dt); %Tridiagonal matrix for advection and diffusion
    
    dSz_inorg = rho_sed*S_resusp.*F_IM.*(-diff([Az; 0])./Vz);  % Dry inorganic particle resuspension source from sediment (kg m-3 day-1)
    Sz = Fi_ad \ (Sz + dSz_inorg);           %Solving new suspended solids profile (advection + diffusion)
    
    dPP = dSz_inorg.*Psz_store;  % PP resuspension source from sediment((kg m-3 day-1)*(mg kg-1) = mg m-3 day-1)
    PPz = Fi_ad \ (PPz + dPP);     %Solving new suspended particulate inorganic P profile (advection + diffusion)
    
    %Chlorophyll, Group 1+2 resuspension (now divided 50/50 between the groups)
    dSz_org = rho_org*S_resusp.*F_LOM.*(-diff([Az; 0])./Vz);  %Dry organic particle resuspension source from sediment (kg m-3 day-1)
    dChl_res = dSz_org.*Chlsz_store;  %Chl a resuspension source from sediment resusp. ((kg m-3 day-1)*(mg kg-1) = mg m-3 day-1);
    
    %Chlorophyll, Group 1
    dChl_growth = Chlz .* R_bioz; %Chl a growth source
    dChl_gross_growth = Chlz .* Growth_bioz; % Gross "biomass" production
    dChl_loss = Chlz .* Loss_bioz;
    dChl = dChl_growth + 0.5*dChl_res; % Total Chl a source (resuspension 50/50 between the two groups, NEW!!!)
    %DOC exudation - evenly in all subpools for now
    DOCz1 = DOCz1 + DOC_excr.*C_per_Chla.*dChl_gross_growth;
    %Phytoplankton respiration
    CO2z = CO2z + PP_resp.*(44.01/12).*C_per_Chla.*dChl_gross_growth;
    O2z = O2z - PP_resp.*(32/12).*C_per_Chla.*dChl_gross_growth;
    %Transfer to POC
    POCz1 = POCz1 + C_per_Chla.*dChl_loss;
    
    Fi_ad = tridiag_HAD_v11([NaN; Kz],w_chl,Vz,Az,dz,dt); %Tridiagonal matrix for advection and diffusion
    Chlz = Fi_ad \ (Chlz + dChl);  %Solving new phytoplankton profile (advection + diffusion) (always larger than background level)
    
    %Chlorophyll, Group 2
    dCz_growth = Cz .* R_bioz_2; %Chl a growth source
    dCz_gross_growth = Cz .* Growth_bioz_2; % Gross "biomass" production
    dCz_loss = Cz .* Loss_bioz_2;
    dCz = dCz_growth + 0.5*dChl_res; % Total Chl a source (resuspension 50/50 between the two groups, NEW!!!)
    %DOC exudation - evenly in all subpools for now
    DOCz1 = DOCz1 + DOC_excr.*C_per_Chla.*dCz_gross_growth;
    %Phytoplankton respiration
    %CO2z = CO2z + PP_resp.*(44.01/12).*C_per_Chla.*Cz;
    CO2z = CO2z + PP_resp.*(44.01/12).*C_per_Chla.*dCz_gross_growth;
    %O2z = O2z - PP_resp.*(32/12).*C_per_Chla.*Cz;
    O2z = O2z - PP_resp.*(32/12).*C_per_Chla.*dCz_gross_growth;
    %Transfer to POC
    POCz1 = POCz1 + C_per_Chla.*dCz_loss;
    
    Fi_ad = tridiag_HAD_v11([NaN; Kz],w_chl_2,Vz,Az,dz,dt); %Tridiagonal matrix for advection and diffusion
    Cz = Fi_ad \ (Cz + dCz);  %Solving new phytoplankton profile (advection + diffusion) (always larger than background level)
    
    %Dissolved inorganic phosphorus
    dP = dDOP - (dChl_gross_growth + dCz_gross_growth)./ Y_cp;
    Pz = Fi \ (Pz + dP); %Solving new dissolved inorganic P profile (diffusion)
    
    dmmuy = 0.5*dSz_inorg;
    dmmuy(1:end-3) = 0;
    Salz = Fi \ (Salz + dmmuy);
    
    %Allochthonous (pools 2 & 3) DOC flocculation
    floc2 = flocc*DOCz2;
    floc3 = flocc*DOCz3;
    DOCz2 = DOCz2 - floc2;
    DOCz3 = DOCz3 - floc3;
    POCz2 = POCz2 + floc2 + floc3 + dPOC_res;
    
    %POC goes to the 2 fastest subpool
    DOCz1 = DOCz1 + dPOC1;
    DOCz2 = DOCz2 + 0.5*dPOC2;
    DOCz3 = DOCz3 + 0.5*dPOC2;
    
    %Particulate organic carbon
    Fi_ad = tridiag_HAD_v11([NaN; Kz],w_chl,Vz,Az,dz,dt);
    
    POCz1 = Fi_ad \ POCz1;  %Solving new POC (advection + diffusion) (always larger than background level)
    POCz2 = Fi_ad \ POCz2;

    %Dissolved organic carbon, Fokema - current version
    DOCz_old = DOCz1+DOCz2+DOCz3;
    Currdate=datevec(tt(i)); %Date
    Date=Currdate(1,2); %Month number
    
    kerroin1 = -log(1-RQuot*(1-exp(-k_DOC1)))/k_DOC1;
    kerroin2 = -log(1-RQuot*(1-exp(-k_DOC2)))/k_DOC2;
    
    %Fokema
    [DOCz_new,DOC1_frac,DOC2_frac,DOC3_frac,~,Daily_BB1,Daily_BB2,Daily_BB3,Daily_PB] = fokema_new(DOCz1,DOCz2,DOCz3,0,Qsw,Tz,0,Date,zz,Qlambda,O2z,kerroin1*k_DOC1,kerroin2*k_DOC2,theta_cold,theta_DOC,theta_DOC_ice);
    DOCz1 = DOC1_frac.*DOCz_new;
    DOCz2 = DOC2_frac.*DOCz_new;
    DOCz3 = DOC3_frac.*DOCz_new;
    
    DOCz1 = Fi \ DOCz1; %Solving new dissolved organic C profile (diffusion)
    DOCz2 = Fi \ DOCz2;
    DOCz3 = Fi \ DOCz3;
    
    %TSA model
    %dDOC = -oc_DOC*qy_DOC*f_par*(1/e_par)*(60*60*24*dt)*Qsw*Attn_z; %photochemical degradation
    %[m2/mg_doc]*[mg_doc/mol_qnt]*[-]*[mol_qnt/J]*[s/day]*[J/s/m2]*[-] = [1/day]
    %DOCz = Fi \ (DOCz + dDOC.*DOCz); %Solving new dissolved organic C profile (diffusion)
    
    %Oxygen production in phytoplankton growth
    kerroin = C_per_Chla*32/12*1.03125;
    dO2_Chl = PQuot*kerroin.*(dChl_gross_growth+dCz_gross_growth);
    
    %DO consumption by DOC degradation
    dO2_DOC = 1/RQuot*32./12*(DOCz_old-DOCz_new);
    O2_old = O2z;
    O2z = max(0,O2z + dO2_Chl - dO2_DOC);
    O2_diff = O2z - O2_old;
    
    Fi_O2 = tridiag_DIF_v11([NaN; Kz_O2],Vz,Az,dz,dt); %Tridiagonal matrix for oxygen diffusion
    O2z = Fi_O2 \ O2z; %Solving new dissolved oxygen profile (diffusion)
    
    %Dissolved inorganic carbon
    
    %Oxygen to CO2 yield factor
    
    if(i==DICpaiva);keyboard;end
    CO2z = max(0,CO2z + 1.375.*(RQuot*dO2_DOC-1/PQuot*dO2_Chl));%-0.75*dO2_Chl
    
    [CO2z,DICz,Hplusz] = DICsystem_bio(DICz,CO2z,Tz,Hplusz,i);
    
    if(i==DICpaiva);keyboard;end
    [CO2z,DICz,Hplusz] = DIC_diffusion2(DICz,Tz,Hplusz,Fi_O2); %Solving new dissolved inorganic carbon profile (diffusion)
    if(i==DICpaiva);keyboard;end
    %----------------------------------------------------------------------
    
    %Sediment-water exchange (DOP source neglected)
    %-porewater to water
    
    PwwFracPtW=ksw*(-diff([Az; 0]))./Vz; %fraction between resuspended porewater and water layer volumes
    EquP1 = (1-PwwFracPtW).*Pz + PwwFracPtW.*Pdz_store; %Mixture of porewater and water
    dPW_up = EquP1-Pz; %"source/sink" for output purposes
    
    %-water to porewater
    PwwFracWtP=ksw./((1-F_sed_sld)*H_sed); %NEW testing 3.8.05; fraction between resuspended (incoming) water and sediment layer volumes
    EquP2 = PwwFracWtP.*Pz + (1-PwwFracWtP).*Pdz_store; %Mixture of porewater and water
    dPW_down = EquP2-Pdz_store; %"source/sink" for output purposes
    
    %-update concentrations
    Pz = EquP1;
    Pdz_store=EquP2;
    
    %Yhdistetään POCit sedimenttiä varten
    POCz = POCz1+POCz2;
    POC1_frac = POCz1./POCz;
    
    %Calculate the thickness ratio of newly settled net sedimentation and mix these
    %two to get new sediment P concentrations in sediment (taking into account particle resuspension)
    delPP_inorg=NaN*ones(Nz,1); %initialize
    delC_inorg=NaN*ones(Nz,1); %initialize
    delC_org=NaN*ones(Nz,1); %initialize
    delC_org2=NaN*ones(Nz,1); %initialize % NEW!!! for chlorophyll group 2
    delC_porg=NaN*ones(Nz,1);
    
    delA=diff([Az; 0]); %Area difference for layer i (OBS: negative)
    meanA=0.5*(Az+[Az(2:end); 0]);
    
    %sedimentation is calculated from "Funnelling-NonFunnelling" difference
    %(corrected 03.10.05)
    delPP_inorg(1)=(0 - PPz(1)*delA(1)./meanA(1))./(dz/(dt*w_s) + 1);
    delC_inorg(1)=(0 - Sz(1)*delA(1)./meanA(1))./(dz/(dt*w_s) + 1);
    delC_org(1)=(0 - Chlz(1)*delA(1)./meanA(1))./(dz/(dt*w_chl) + 1);
    delC_org2(1)= (0 - Cz(1)*delA(1)./meanA(1))./(dz/(dt*w_chl_2) + 1); % NEW!!!  for chlorophyll group 2
    delC_porg(1)= (0 - POCz(1)*delA(1)./meanA(1))./(dz/(dt*w_chl) + 1);
    
    for ii=2:Nz
        delPP_inorg(ii)=(delPP_inorg(ii-1) - PPz(ii)*delA(ii)./meanA(ii))./(dz/(dt*w_s) + 1); %(mg m-3)
        delC_inorg(ii)=(delC_inorg(ii-1) - Sz(ii)*delA(ii)./meanA(ii))./(dz/(dt*w_s) + 1); %(kg m-3)
        delC_org(ii)=(delC_org(ii-1) - Chlz(ii)*delA(ii)./meanA(ii))./(dz/(dt*w_chl) + 1); %(mg m-3)
        delC_org2(ii)=(delC_org2(ii-1) - Cz(ii)*delA(ii)./meanA(ii))./(dz/(dt*w_chl_2) + 1); %(mg m-3) % NEW!!! for chlorophyll group 2
        delC_porg(ii)=(delC_porg(ii-1) - POCz(ii)*delA(ii)./meanA(ii))./(dz/(dt*w_chl) + 1); %(mg m-3)
    end
    
    H_netsed_inorg=max(0, (Vz./(-diff([Az; 0]))).*delC_inorg./rho_sed - F_IM.*S_resusp); %inorganic(m day-1, dry), always positive
    H_netsed_org=max(0, (Vz./(-diff([Az; 0]))).*(delC_org+delC_org2)./(F_OM*Y_cp*rho_org) - F_LOM.*S_resusp);
    %organic (m day-1, dry), always positive,  NEW!!! for chlorophyll group 2
    H_netsed_porg=max(0, (Vz./(-diff([Az; 0]))).*(delC_porg./(F_OC*rho_org)) - F_NLOM.*S_resusp); %PK
    
    %H_totsed=H_netsed_org + H_netsed_inorg;  %total (m day-1), always positive
    H_totsed=H_netsed_org + H_netsed_inorg + H_netsed_porg; %total (m day-1), always positive
    
    %F_IM_NewSed=F_IM;
    F_LOM_NewSed=F_LOM;
    F_NLOM_NewSed=F_NLOM;
    inx=find(H_totsed>0);
    %F_IM_NewSed(inx)=H_netsed_inorg(inx)./H_totsed(inx); %volume fraction of inorganic matter in net settling sediment
    F_LOM_NewSed(inx)=H_netsed_org(inx)./H_totsed(inx);
    F_NLOM_NewSed(inx)=H_netsed_porg(inx)./H_totsed(inx);
    
    NewSedFrac = min(1, H_totsed./(F_sed_sld*H_sed)); %Fraction of newly fallen net sediment of total active sediment depth, never above 1
    NewSedFrac_inorg = min(1, H_netsed_inorg./(F_IM.*F_sed_sld*H_sed)); %Fraction of newly fallen net inorganic sediment of total active sediment depth, never above 1
    NewSedFrac_org = min(1, H_netsed_org./(F_LOM.*F_sed_sld*H_sed)); %Fraction of newly fallen net organic sediment of total active sediment depth, never above 1
    NewSedFrac_porg = min(1, H_netsed_porg./(F_NLOM.*F_sed_sld*H_sed)); %Fraction of newly fallen net nonliving organic sed. of total active sed. depth, never above 1
    
    %Psz_store: %P conc. in inorganic sediment particles (mg kg-1 dry w.)
    Psz_store = (1-NewSedFrac_inorg).*Psz_store + NewSedFrac_inorg.*PPz./Sz; %(mg kg-1)
    
    %Update counters
    Sedimentation_counter = Sedimentation_counter + Vz.*(delC_inorg + (delC_org+delC_org2)./(F_OM*Y_cp)); %Inorg.+Org. (kg)
    POC_sedim_counter = POC_sedim_counter + Vz.*delC_porg./F_OC; %(kg nonliving organic matter) %lasketaan siis epäelävänä orgaanisena _aineena_
    Resuspension_counter = Resuspension_counter + Vz.*(dSz_inorg + dSz_org); %Inorg.+Org. (kg)
    POC_resusp_counter = POC_resusp_counter + Vz.*dSz_porg; %Org. (kg)
    
    %Chlsz_store (for group 1+2): %Chl a conc. in sediment particles (mg kg-1 dry w.)
    Chlsz_store = (1-NewSedFrac_org).*Chlsz_store + NewSedFrac_org.*F_OM*Y_cp; %(mg kg-1)
    %Subtract degradation to P in pore water
    
    O2_fact = O2z./(500+O2z);
    
    %Ei lasketa suoraan hajoamista, vaan otetaan happiehto:
    Chlz_seddeg_max = 0.85*k_twty * Y_cp*F_OM*rho_org* F_LOM.*O2_fact.* theta_m.^(Tz-20); %(mg Chl)(m-3 dry sed)
    %POCin suora hajotus
    POCz_seddeg_max = k_POC_twty*F_OC*rho_org * F_NLOM.*O2_fact.*POC_sed_t_fctr; %(mg C)(m-3 dry sed)
    
    POP_seddeg_max = 0.11765*k_POC_twty*POPz_sed.*O2_fact.*POC_sed_t_fctr; %(mg POP)(m-3 dry sed)
    %Ennen POCz_seddeg_max = 0.85*1.7*k_POC_twty ja POP_seddeg_max =
    %0.17*k_POC_twty. Nyt etukertoimet mukana k_POC_twty:ssä, joten POP-kerroin pitää skaalata.
    
    V_dry_sed = H_sed*(-diff([Az; 0]))*F_sed_sld;
    
    %mass of max. degrading Chl
    m_Chlz_seddeg_max = Chlz_seddeg_max.*V_dry_sed; %mg Chl
    %mass of total oxygen in the corresponding water layer
    m_O2z = O2z.*Vz; %mg O2
    %mass of consumed oxygen in maximum case
    m_dO2_chl_sed_max = 1/RQuot*110*m_Chlz_seddeg_max; %mg O2
    dO2_chl_sed = NaN*ones(Nz,1); %(mg O2)(m-3 water)
    for h = 1:length(zz);
        if(m_O2z(h) < m_dO2_chl_sed_max(h))
            dO2_chl_sed(h) = O2z(h);
        else
            dO2_chl_sed(h) = m_dO2_chl_sed_max(h)./Vz(h);
        end
    end
    Chlz_seddeg = dO2_chl_sed.*Vz./(1/RQuot*110*V_dry_sed);
    O2z = O2z-dO2_chl_sed;
    
    m_POCz_seddeg_max = POCz_seddeg_max.*V_dry_sed; %mg C
    m_O2z = O2z.*Vz; %mg O2
    m_dO2_poc_sed_max = 1/RQuot*32/12*m_POCz_seddeg_max; %mg O2
    dO2_sed = NaN*ones(Nz,1);
    for h = 1:length(zz);
        if(m_O2z(h) < m_dO2_poc_sed_max(h))
            dO2_sed(h) = O2z(h);
        else
            dO2_sed(h) = m_dO2_poc_sed_max(h)./Vz(h);
        end
    end
    POCz_seddeg = dO2_sed.*Vz./(1/RQuot*32/12*V_dry_sed); % (mg C) (m-3 dry organic sediment)
    
    POPz_seddeg = POCz_seddeg./POCz_seddeg_max.*POP_seddeg_max;
    
    %------
    O2z = O2z-dO2_sed;
    
    dO2_SOD = dO2_chl_sed+dO2_sed;
    
    CO2z = max(0,CO2z + RQuot*1.375.*dO2_SOD);
    
    [CO2z,DICz,Hplusz] = DICsystem_bio(DICz,CO2z,Tz,Hplusz,i);
    
    if(i==DICpaiva);keyboard;end
    %calculate new VOLUME fraction of inorganic particles of total dry sediment
    F_LOM = min(1,F_LOM.*(1 - Chlz_seddeg/(Y_cp*F_OM*rho_org)) ).*(1-NewSedFrac) + F_LOM_NewSed.*NewSedFrac;
    F_NLOM = min(1, F_NLOM.*(1 - POCz_seddeg/(F_OC*rho_org)) ).*(1-NewSedFrac) + F_NLOM_NewSed.*NewSedFrac;
    F_IM = 1 - F_LOM- F_NLOM;
    
    P_frac_new = 1./(POC1_frac./(C_per_P_pp) + (1-POC1_frac)./C_per_P_alloc);
    POPz_sed = POPz_sed - POPz_seddeg;
    %mdd = POPz_sed;
    C_P_org_sed_old = (F_OC*rho_org * F_NLOM)./POPz_sed;
    
    %Pdz_store=Pdz_store + Chlz_seddeg .* F_sed_sld/((1-F_sed_sld)*Y_cp) + POCz_seddeg .* F_sed_sld./((1-F_sed_sld).*C_per_P_alloc);%P_frac_new%C_P_org_sed);
    %C_P_org_sed = (1-NewSedFrac_porg).*C_P_org_sed + NewSedFrac_porg.*P_frac_new;
    Pdz_store=Pdz_store + Chlz_seddeg .* F_sed_sld/((1-F_sed_sld)*Y_cp) + POPz_seddeg .* F_sed_sld./((1-F_sed_sld));
    C_P_org_sed = 1./((1-NewSedFrac_porg)./C_P_org_sed_old + NewSedFrac_porg./P_frac_new);
    POPz_sed = (F_OC*rho_org * F_NLOM)./C_P_org_sed;
    
    %== P-partitioning in sediments==
    VolFrac=1./(1+(1-F_sed_sld)./(F_sed_sld*F_IM)); %volume fraction: inorg sed. / (inorg.sed + pore water)
    TIP_sed =rho_sed*VolFrac.*Psz_store + (1-VolFrac).*Pdz_store; %total inorganic P in sediments (mg m-3)
    
    [Pdz_store, Psz_store]=Ppart(VolFrac,TIP_sed,Psat_L,Fmax_L_sed,rho_sed,Fstable);
    
    
    %Oxygen surface flux
    if(IceIndicator==0)
        [O2z(1),O2flux,O2_eq,K0_O2] = oxygenflux(O2z(1),Wt(i,6),Wt(i,5),Tz(1),dz);
    else
        O2flux = 0;
        O2_eq = NaN;
        K0_O2 = NaN;
    end
    
    
    %Carbon dioxide surface flux
    if(IceIndicator==0)
        [CO2z,surfflux,CO2_eq,K0,CO2_ppm] = carbondioxideflux_new(CO2z,Wt(i,6),Wt(i,5),Tz(1),dz,tt(i),Az,Vz,CO2air);
        if(CO2z(1)<0);disp(num2str(i));disp(num2str(CO2z(1)));keyboard;end
        
        [CO2z,DICz,Hplusz] = DICsystem_new(DICz,CO2z,Tz,Tz,Hplusz,i);
        
    else
        surfflux=0;
        CO2_eq = NaN;
        K0 = NaN;
        CO2_ppm = NaN;
    end
    if(i==DICpaiva);keyboard;end
    % Inflow calculation
    % Inflw(:,1) Inflow volume (m3 day-1)
    % Inflw(:,2) Inflow temperature (deg C)
    % Inflw(:,3) Inflow chlorophyll (group 2) concentration (-)
    % Inflw(:,4) Inflow sedimenting tracer (or suspended inorganic matter) concentration (kg m-3)
    % Inflw(:,5) Inflow total phosphorus (TP) concentration  (incl. DOP & Chla) (mg m-3)
    % Inflw(:,6) Inflow dissolved organic phosphorus (DOP) concentration (mg m-3)
    % Inflw(:,7) Inflow chlorophyll (group 1) concentration (mg m-3)
    % Inflw(:,8) Inflow DOC concentration (mg m-3)
    % Inflw(:,9) Inflow DIC concentration (mg m-3)
    % Inflw(:,10) Inflow O2 concentration (mg m-3)
    % Inflw(:,11) Inflow POC concentration (mg m-3)
    
    if (river_inflow_switch==1)
        Iflw = I_scV * Inflw(i,1); % (scaled) inflow rate
        Iflw_T = I_scT + Inflw(i,2); %(adjusted) inflow temperature
        if (Iflw_T<Tf) %negative temperatures changed to Tf
            Iflw_T=Tf;
        end
        Iflw_C = I_scC * Inflw(i,3); %(scaled) inflow C concentration
        Iflw_S = I_scS * Inflw(i,4); %(scaled) inflow S concentration
        Iflw_TP = I_scTP * Inflw(i,5); %(scaled) inflow TP concentration (incl. DOP & Chla)
        Iflw_DOP = I_scDOP * Inflw(i,6); %(scaled) inflow DOP concentration
        Iflw_Chl = I_scChl * Inflw(i,7); %(scaled) inflow Chl a concentration
        Iflw_DOC = I_scDOC * Inflw(i,8); %(scaled) inflow DOC concentration
        Iflw_DIC = I_scDIC * Inflw(i,9); %(scaled) inflow DIC concentration
        Iflw_O2 = I_scDO * Inflw(i,10); %(scaled) inflow O2 concentration
        Iflw_POC = I_scDOC * Inflw(i,11); %(scaled) inflow POC concentration
        Inflw_Hplus = Inflw(i,12);
        Iflw_Sal = 0 * Inflw(i,11);
        
        %Added suspended solids correction: minimum allowed P bioavailability factor is 0.1
        if any((1-(Iflw_DOP+(Iflw_Chl+Iflw_C)./Y_cp+Iflw_POC./C_per_P_alloc)./Iflw_TP-(Iflw_S*Fstable)./Iflw_TP)<0.1); % NEW!!!!
            Iflw_S_dum = (1 - (Iflw_DOP+(Iflw_Chl+Iflw_C)./Y_cp+Iflw_POC./C_per_P_alloc)./Iflw_TP - 0.1).*(Iflw_TP./Fstable); %NEW!!!
            SS_decr=SS_decr+(Iflw_S-Iflw_S_dum)*Iflw;
            Iflw_S=Iflw_S_dum;%disp(num2str(i));
        end
        %keyboard
        if any((Iflw_TP-Iflw_DOP-(Iflw_Chl+Iflw_C)./Y_cp-Iflw_POC./C_per_P_alloc-Iflw_S*Fstable)<0)  %NEW!!!
            error('Sum of DOP, inactive PP, and P contained in Chl a (both groups)and POC in inflow cannot be larger than TP')
        end
        %keyboard
        if(Iflw>0)
            if (isnan(Iflw_T))
                lvlD=0;
                Iflw_T=Tz(1);
            else
                if(density_switch == 0)
                    rho = polyval(ies80,max(0,Tz(:)))+min(Tz(:),0);	% Density (kg/m3)
                else
                    rho = sw_dens(Salz,max(0,Tz(:)),0) + min(Tz(:),0);
                end
                rho_Iflw=polyval(ies80,max(0,Iflw_T))+min(Iflw_T,0);
                lvlG=find(rho>=rho_Iflw);
                if (isempty(lvlG))
                    lvlG=length(rho);
                end
                lvlD=min(zz(floor(end/2)),zz(lvlG(1))); %level zz above which inflow is put
                lvlD_DIC = zz(end);
            end %if isnan...
            
            Tz_IC_old = Tz;
            %Changes in properties due to inflow
            dummy=IOflow_v11(dz, zz, Vz, Tz, lvlD, Iflw, Iflw_T);
            Tz=dummy; %Temperature
            
            [CO2z,DICz,Hplusz] = DICsystem_new(DICz,CO2z,Tz_IC_old,Tz,Hplusz,i);
            if(i==DICpaiva);keyboard;end
            dummy=IOflow_v11(dz, zz, Vz, Sz, lvlD, Iflw, Iflw_S);
            Sz=dummy; %Sedimenting tracer
            
            dummy=IOflow_v11(dz, zz, Vz, DOPz, lvlD, Iflw, Iflw_DOP);
            DOPz=dummy; %Particulate organic P
            
            TIPz=Pz + PPz; % Total inorg. phosphorus (excl. Chla and DOP) in the water column (mg m-3)
            dummy=IOflow_v11(dz, zz, Vz, TIPz, lvlD, Iflw, Iflw_TP-((Iflw_Chl+Iflw_C)./Y_cp)-Iflw_DOP); %NEW!!!
            TIPz=dummy; %Total inorg. phosphorus (excl. Chla and DOP)
            
            
            %== P-partitioning in water==
            [Pz, ~]=Ppart(Sz./rho_sed,TIPz,Psat_L,Fmax_L,rho_sed,Fstable);
            PPz=TIPz-Pz;
            
            dummy=IOflow_v11(dz, zz, Vz, Cz, lvlD, Iflw, Iflw_C);
            Cz=dummy; %Chlorophyll (group 2)
            
            dummy=IOflow_v11(dz, zz, Vz, Chlz, lvlD, Iflw, Iflw_Chl);
            Chlz=dummy; %Chlorophyll (group 1)
            
            %DOC inflow to subpools 2&3
            dmmfctr = 0.1875;
            
            dummy=IOflow_v11(dz, zz, Vz, DOCz2, lvlD, Iflw, dmmfctr*Iflw_DOC);
            DOCz2=dummy; %DOC
            
            dummy=IOflow_v11(dz, zz, Vz, DOCz3, lvlD, Iflw, (1-dmmfctr)*Iflw_DOC);
            DOCz3=dummy; %DOC
            
            dummy=IOflow_v11(dz, zz, Vz, O2z, lvlD, Iflw, Iflw_O2);
            O2z=dummy; %O2
            
            dummy=IOflow_v11(dz, zz, Vz, POCz1, lvlD, Iflw, 0);
            POCz1 = dummy;%POC
            
            dummy=IOflow_v11(dz, zz, Vz, POCz2, lvlD, Iflw, Iflw_POC);
            POCz2 = dummy;
            
            dummy=IOflow_v11(dz, zz, Vz, Salz, lvlD, Iflw, Iflw_Sal);
            Salz=dummy; %Salinity
            
            %[CO2z,DICz,Hplusz] = IOflow_DIC(dz, zz, Vz, Tz, DICz, Hplusz, lvlD, Iflw, Tz(1), DICz(1), Hplusz(1));
            [CO2z,DICz,Hplusz] = IOflow_DIC(dz, zz, Vz, Tz, DICz, Hplusz, lvlD, Iflw, Iflw_T, Iflw_DIC, Inflw_Hplus);
        else
            lvlD=NaN;
        end %if(Iflw>0)
        
    else
        %These are needed only in MixStat counting (P budget)
        Iflw=0; % (scaled) inflow rate
        %Iflw_T = NaN; %(adjusted) inflow temperature
        Iflw_C = NaN; %(scaled) inflow C concentration
        Iflw_S = NaN; %(scaled) inflow S concentration
        Iflw_TP = NaN; %(scaled) inflow TP concentration (incl. DOP & Chla)
        Iflw_DOP = NaN; %(scaled) inflow DOP concentration
        Iflw_Chl = NaN; %(scaled) inflow Chl a concentration
        %Iflw_DOC = NaN; %(scaled) inflow DOC concentration
        %Iflw_DIC = NaN; %(scaled) inflow DIC concentration
        %Iflw_O2 = NaN; %(scaled) inflow O2 concentration
        %Iflw_POC = NaN; %(scaled) inflow POC concentration
        lvlD=NaN;
    end  %if (river_inflow_switch==1)
    
    if(i==DICpaiva);keyboard;end
    % Convective mixing adjustment (mix successive layers until stable density profile,  no summer/autumn turnover here!)
    
    [Tz,Cz,Sz,Pz,Chlz,PPz,DOPz,DOCz1,DOCz2,DOCz3,DICz,O2z,POCz1,POCz2,Hplusz,CO2z] = convection_v12_1aIC(Tz,Cz,Sz,Pz,Chlz,PPz,DOPz,DOCz1,DOCz2,DOCz3,DICz,O2z,POCz1,POCz2,Hplusz,CO2z,Tprof_prev,Vz,Cw,f_par,lambdaz_wtot_avg,zz,swa_b0,tracer_switch,0,i);
    if(i==naytto);display(['Inflow ',num2str(Tz')]);end
    
    if(i==DICpaiva);display('Juuri ennen wind mixingiä');keyboard;end
    if (IceIndicator==0)
        if(density_switch == 0)
            rho = polyval(ies80,max(0,0.5*(T_testi(:)+Tz(:))))+min(0.5*(T_testi(:)+Tz(:)),0);
        else
            rho = sw_dens(Salz,max(0,Tz(:)),0) + min(Tz(:),0);
        end
        TKE=C_shelter*Az(1)*sqrt(tau^3/rho(1))*(24*60*60*dt); %Turbulent kinetic energy (J day-1) over the whole lake
        Tz_wmIC_old = Tz;
        %Wind mixing
        WmixIndicator=1;
        Bef_wind=sum(diff(rho)==0); %just a watch variable
        while (WmixIndicator==1)
            d_rho=diff(rho);
            inx=find(d_rho>0);
            if (isempty(inx)==0); %if water column not already fully mixed
                zb=inx(1);
                MLD=dz*zb; %mixed layer depth
                dD=d_rho(zb); %density difference
                Zg=sum( Az(1:zb+1) .* zz(1:zb+1) ) / sum(Az(1:zb+1)); %Depth of center of mass of mixed layer
                V_weight=Vz(zb+1)*sum(Vz(1:zb))/(Vz(zb+1)+sum(Vz(1:zb)));
                POE=(dD*g*V_weight*(MLD + dz/2 - Zg));
                KP_ratio=TKE/POE;
                if (KP_ratio>=1)
                    Tmix=sum( Vz(1:zb+1).*Tz(1:zb+1) ) / sum(Vz(1:zb+1));
                    Tz(1:zb+1)=Tmix;
                    
                    Cmix=sum( Vz(1:zb+1).*Cz(1:zb+1) ) / sum(Vz(1:zb+1));
                    Cz(1:zb+1)=Cmix;
                    
                    Smix=sum( Vz(1:zb+1).*Sz(1:zb+1) ) / sum(Vz(1:zb+1));
                    Sz(1:zb+1)=Smix;
                    
                    Pmix=sum( Vz(1:zb+1).*Pz(1:zb+1) ) / sum(Vz(1:zb+1));
                    Pz(1:zb+1)=Pmix;
                    
                    Chlmix=sum( Vz(1:zb+1).*Chlz(1:zb+1) ) / sum(Vz(1:zb+1));
                    Chlz(1:zb+1)=Chlmix;
                    
                    PPmix=sum( Vz(1:zb+1).*PPz(1:zb+1) ) / sum(Vz(1:zb+1));
                    PPz(1:zb+1)=PPmix;
                    
                    DOPmix=sum( Vz(1:zb+1).*DOPz(1:zb+1) ) / sum(Vz(1:zb+1));
                    DOPz(1:zb+1)=DOPmix;
                    
                    DOC1mix=sum( Vz(1:zb+1).*DOCz1(1:zb+1) ) / sum(Vz(1:zb+1));
                    DOCz1(1:zb+1)=DOC1mix;
                    
                    DOC2mix=sum( Vz(1:zb+1).*DOCz2(1:zb+1) ) / sum(Vz(1:zb+1));
                    DOCz2(1:zb+1)=DOC2mix;
                    
                    DOC3mix=sum( Vz(1:zb+1).*DOCz3(1:zb+1) ) / sum(Vz(1:zb+1));
                    DOCz3(1:zb+1)=DOC3mix;
                    
                    O2mix=sum( Vz(1:zb+1).*O2z(1:zb+1) ) / sum(Vz(1:zb+1));
                    O2z(1:zb+1)=O2mix;
                    
                    POC1mix=sum( Vz(1:zb+1).*POCz1(1:zb+1) ) / sum(Vz(1:zb+1));
                    POCz1(1:zb+1)=POC1mix;
                    
                    POC2mix=sum( Vz(1:zb+1).*POCz2(1:zb+1) ) / sum(Vz(1:zb+1));
                    POCz2(1:zb+1)=POC2mix;
                    
                    Salmix=sum( Vz(1:zb+1).*Salz(1:zb+1) ) / sum(Vz(1:zb+1));
                    Salz(1:zb+1)=Salmix;
                    
                    if(density_switch == 0)
                        rho = polyval(ies80,max(0,Tz(:))) + min(Tz(:),0);
                    else
                        rho = sw_dens(Salz,max(0,Tz(:)),0) + min(Tz(:),0);
                    end
                    TKE=TKE-POE;
                else %if KP_ratio < 1, then mix with the remaining TKE part of the underlying layer
                    Tmix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*Tz(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
                    Tz(1:zb)=Tmix;
                    Tz(zb+1)=KP_ratio*Tmix + (1-KP_ratio)*Tz(zb+1);
                    
                    Cmix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*Cz(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
                    Cz(1:zb)=Cmix;
                    Cz(zb+1)=KP_ratio*Cmix + (1-KP_ratio)*Cz(zb+1);
                    
                    Smix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*Sz(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
                    Sz(1:zb)=Smix;
                    Sz(zb+1)=KP_ratio*Smix + (1-KP_ratio)*Sz(zb+1);
                    
                    Pmix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*Pz(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
                    Pz(1:zb)=Pmix;
                    Pz(zb+1)=KP_ratio*Pmix + (1-KP_ratio)*Pz(zb+1);
                    
                    Chlmix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*Chlz(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
                    Chlz(1:zb)=Chlmix;
                    Chlz(zb+1)=KP_ratio*Chlmix + (1-KP_ratio)*Chlz(zb+1);
                    
                    PPmix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*PPz(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
                    PPz(1:zb)=PPmix;
                    PPz(zb+1)=KP_ratio*PPmix + (1-KP_ratio)*PPz(zb+1);
                    
                    DOPmix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*DOPz(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
                    DOPz(1:zb)=DOPmix;
                    DOPz(zb+1)=KP_ratio*DOPmix + (1-KP_ratio)*DOPz(zb+1);
                    
                    DOC1mix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*DOCz1(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
                    DOCz1(1:zb)=DOC1mix;
                    DOCz1(zb+1)=KP_ratio*DOC1mix + (1-KP_ratio)*DOCz1(zb+1);
                    
                    DOC2mix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*DOCz2(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
                    DOCz2(1:zb)=DOC2mix;
                    DOCz2(zb+1)=KP_ratio*DOC2mix + (1-KP_ratio)*DOCz2(zb+1);
                    
                    DOC3mix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*DOCz3(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
                    DOCz3(1:zb)=DOC3mix;
                    DOCz3(zb+1)=KP_ratio*DOC3mix + (1-KP_ratio)*DOCz3(zb+1);
                    
                    O2mix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*O2z(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
                    O2z(1:zb)=O2mix;
                    O2z(zb+1)=KP_ratio*O2mix + (1-KP_ratio)*O2z(zb+1);
                    
                    POC1mix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*POCz1(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
                    POCz1(1:zb)=POC1mix;
                    POCz1(zb+1)=KP_ratio*POC1mix + (1-KP_ratio)*POCz1(zb+1);
                    
                    POC2mix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*POCz2(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
                    POCz2(1:zb)=POC2mix;
                    POCz2(zb+1)=KP_ratio*POC2mix + (1-KP_ratio)*POCz2(zb+1);
                    
                    Salmix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*Salz(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
                    Salz(1:zb)=Salmix;
                    Salz(zb+1)=KP_ratio*Salmix + (1-KP_ratio)*Salz(zb+1);
                    
                    if(density_switch == 0)
                        %rho = polyval(ies80,max(0,Tz(:))) + min(Tz(:),0);
                        rho = polyval(ies80,max(0,0.5*(T_testi(:)+Tz(:))))+min(0.5*(T_testi(:)+Tz(:)),0);
                    else
                        rho = sw_dens(Salz,max(0,Tz(:)),0) + min(Tz(:),0);
                    end
                    TKE=0;
                    WmixIndicator=0;
                    % DIC system mixing calculated separately for the whole mixing layer
                    [CO2z,DICz,Hplusz] = DIC_windmix(DICz,Tz_wmIC_old,Tz,Hplusz,Vz,zb,KP_ratio);
                    %-----------------------------
                end %if (KP_ratio>=1)
            else
                WmixIndicator=0;
                %----------------------------------
                [CO2z,DICz,Hplusz] = DIC_windmix(DICz,Tz_wmIC_old,Tz,Hplusz,Vz,zb,KP_ratio);
                %----------------------------------
            end %if water column (not) already mixed
        end %while
        if(i==DICpaiva);display(['Wind mixingin jälkeen; KP_ratio on ', num2str(KP_ratio)]);keyboard;end
        Aft_wind=sum(diff(rho)==0); %just a watch variable
        Tz_IC_old = Tz;
    else % ice cover module
        Tz_IC_old = Tz;
        
        XE_surf=(Tz(1)-Tf) * Cw * dz; %Daily heat accumulation into the first water layer (J m-2)
        Tz(1)=Tf; %Ensure that temperature of the first water layer is kept at freezing point
        TKE=0; %No energy for wind mixing under ice
        
        if (Wt(i,3)<Tf) %if air temperature is below freezing
            %Calculate ice surface temperature (Tice)
            if(WEQs==0) %if no snow
                alfa=1/(10*Hi);
                dHsi=0;
            else
                K_snow=2.22362*(rho_snow/1000)^1.885; %Yen (1981)
                alfa=(K_ice/K_snow)*(((rho_fw/rho_snow)*WEQs)/Hi);
                %Slush/snow ice formation (directly to ice)
                dHsi=max([0, Hi*(rho_ice/rho_fw-1)+WEQs]);
                Hsi=Hsi+dHsi;
            end
            Tice=(alfa*Tf+Wt(i,3))/(1+alfa);
            %if(i==370);keyboard;end
            %Ice growth by Stefan's law
            Hi_new=sqrt((Hi+dHsi)^2+(2*K_ice/(rho_ice*L_ice))*(24*60*60)*(Tf-Tice));
            %snow fall
            dWEQnews=0.001*Wt(i,7); %mm->m
            if(Wt(i,3)>-1 && Wt(i,7) > 0.15);dWEQnews = 0;end
            dWEQs=dWEQnews-dHsi*(rho_ice/rho_fw); % new precipitation minus snow-to-snowice in snow water equivalent
            %if(i>63);keyboard;end
            dHsi=0; %reset new snow ice formation
        else %if air temperature is NOT below freezing
            Tice=Tf;    %ice surface at freezing point
            dWEQnews=0; %No new snow
            if (WEQs>0)
                %snow melting in water equivalents
                %dWEQs=-max([0, (60*60*24)*(((1-IceSnowAttCoeff)*Qsw)+Qlw+Qsl)/(rho_fw*L_ice)]);
                dWEQs=-max([0, (60*60*24)*(SnowAttCoeff*Qsw+Qlw+Qsl)/(rho_fw*L_ice)]);
                %if(i==460);keyboard;end
                if ((WEQs+dWEQs)<0) %if more than all snow melts...
                    Hi_new=Hi+(WEQs+dWEQs)*(rho_fw/rho_ice); %...take the excess melting from ice thickness
                else
                    Hi_new=Hi; %ice does not melt until snow is melted away
                end
                %keyboard
            else
                %total ice melting
                dWEQs=0;
                %Hi_new=Hi-max([0, (60*60*24)*(((1-IceSnowAttCoeff)*Qsw)+Qlw+Qsl)/(rho_ice*L_ice)]);
                Hi_new=Hi-max([0, (60*60*24)*(IceSnowAttCoeff*Qsw+Qlw+Qsl)/(rho_ice*L_ice)]);
                %snow ice part melting
                %Hsi=Hsi-max([0, (60*60*24)*(((1-IceSnowAttCoeff)*Qsw)+Qlw+Qsl)/(rho_ice*L_ice)]);
                Hsi=Hsi-max([0, (60*60*24)*(IceSnowAttCoeff*Qsw+Qlw+Qsl)/(rho_ice*L_ice)]);
                if (Hsi<=0)
                    Hsi=0;
                end
            end %if there is snow or not
        end %if air temperature is or isn't below freezing
        
        %Update ice and snow thicknesses
        Hi=Hi_new-(XE_surf/(rho_ice*L_ice)); %new ice thickness (minus melting due to heat flux from water)
        XE_surf=0; %reset energy flux from water to ice (J m-2 per day)
        WEQs=WEQs+dWEQs; %new snow water equivalent
        %keyboard
        if(Hi<Hsi)
            Hsi=max(0,Hi);    %to ensure that snow ice thickness does not exceed ice thickness
            %(if e.g. much ice melting much from bottom)
        end
        
        rho_new_snow = 100/(119.17*0.5)*(67.92+51.25*exp(Wt(i,3)/2.59));
        
        if(WEQs<=0)
            WEQs=0; %excess melt energy already transferred to ice above
            rho_snow=250;%rho_new_snow;
        else
            %Update snow density as weighed average of old and new snow densities
            rho_snow=rho_snow*(WEQs-dWEQnews)/WEQs + rho_new_snow*dWEQnews/WEQs;
            if (snow_compaction_switch==1)
                %snow compaction
                if (Wt(i,3)<Tf) %if air temperature is below freezing
                    rhos=1e-3*rho_snow; %from kg/m3 to g/cm3
                    delta_rhos=24*rhos*C1*(0.5*WEQs)*exp(-C2*rhos)*exp(-0.08*(Tf-0.5*(Tice+Wt(i,3))));
                    rho_snow=min([rho_snow+1e+3*delta_rhos, max_rho_snow]);  %from g/cm3 back to kg/m3
                else
                    if(Wt(i,3)>0.5 && Wt(i,7)>0.25)
                        rho_snow=max_rho_snow;
                    else
                        rho_snow_prev = rho_snow;
                        rho_snow = min(max_rho_snow,100/(119.17*0.45)*(67.92+51.25*exp(Wt(i,3)/2.59)));%rho_snow=350;
                        rho_snow = max(rho_snow_prev,rho_snow);
                    end
                end
            end
        end
        
        if(Hi<=0)
            IceIndicator=0;
            disp(['Ice-off, ' datestr(datenum(M_start)+i-1)])
            XE_melt=(-Hi-(WEQs*rho_fw/rho_ice))*rho_ice*L_ice/(24*60*60);
            %(W m-2) snow part is in case ice has melted from bottom leaving some snow on top (reducing XE_melt)
            Hi=0;
            WEQs=0;
            Tice=NaN;
            DoM(pp)=i;
            pp=pp+1;
        end
    end %of ice cover module
    
    
    %== P-partitioning in water==
    TIPz=Pz + PPz; % Total inorg. phosphorus (excl. Chla and DOP) in the water column (mg m-3)
    [Pz, ~]=Ppart(Sz./rho_sed,TIPz,Psat_L,Fmax_L,rho_sed,Fstable);
    PPz=TIPz-Pz;
    
    % Dissolved oxygen saturation
    [O2_sat_rel, O2_sat_abs] = relative_oxygen(O2z,Tz,Wt(i,5),dz);
    
    %Initial freezing
    Supercooled=find(Tz<Tf);
    if (isempty(Supercooled)==0)
        %===NEW!!! (040707)
        if(Supercooled(1)~=1); disp('NOTE: non-surface subsurface supercooling'); end;
        InitIceEnergy=sum((Tf-Tz(Supercooled)).*Vz(Supercooled)*Cw);
        HFrazil=HFrazil+(InitIceEnergy/(rho_ice*L_ice))/Az(1);
        Tz(Supercooled)=Tf;
        
        if ((IceIndicator==0)&&(HFrazil > Frazil2Ice_tresh))
            IceIndicator=1;
            Hi=Hi+HFrazil;
            HFrazil=0;
            DoF(qq)=i;
            disp(['Ice-on, ' datestr(datenum(M_start)+i-1)])
            qq=qq+1;
        end
        
        if (IceIndicator==1)
            Hi=Hi+HFrazil;
            HFrazil=0;
        end
        Tz(1)=Tf; %set temperature of the first layer to freezing point
        %======================
        
    end
    
    %DIC-partitioning in water
    [CO2z,DICz,Hplusz] = DICsystem_new(DICz,CO2z,Tz_IC_old,Tz,Hplusz,i);
    
    if(i==DICpaiva);keyboard;end
    
    %[CO2z,DICz,Hplusz] = pHinit(DIC0,CO20,T0,Hplus0,DICz,CO2z,Tz,Hplusz);
    
    % Calculate pycnocline depth
    pycno_thres=0.1;  %treshold density gradient value (kg m-3 m-1)
    if(density_switch == 0)
        rho = polyval(ies80,max(0,Tz(:))) + min(Tz(:),0);
    else
        rho = sw_dens(Salz,max(0,Tz(:)),0) + min(Tz(:),0);
    end
    dRdz = [NaN; abs(diff(rho))];
    di=find((dRdz<(pycno_thres*dz)) | isnan(dRdz));
    %dRdz(di)=NaN;
    %TCz = nansum(zz .* dRdz) ./ nansum(dRdz);
    dRdz(di)=0; %modified for MATLAB version 7
    TCz = sum(zz .* dRdz) ./ sum(dRdz);
    
    %vector with S_res_epi above, and S_res_hypo below the pycnocline
    inx=find(zz <= TCz);
    S_resusp(inx)=S_res_epi;
    inx=find(zz > TCz);
    S_resusp(inx)=S_res_hypo;
    
    if (IceIndicator==1)
        S_resusp(:)=S_res_hypo;  %only hypolimnetic type of resuspension allowed under ice
    end
    
    if( isnan(TCz) && (IceIndicator==0) )
        S_resusp(:)=S_res_epi;   %if no pycnocline and open water, resuspension allowed from top to bottom
    end
    if(i==naytto);disp(num2str(Tz'));end
    
    DOCz = DOCz1 + DOCz2 + DOCz3;
    
    POCz = POCz1+POCz2;
    
    % Output matrices
    Qst(:,i) = [Qsw Qlw Qsl]';
    Kzt(:,i) = [0;Kz];
    Tzt(:,i) = Tz;
    Czt(:,i) = Cz;
    Szt(:,i) = Sz;
    Pzt(:,i) = Pz;
    Chlzt(:,i) = Chlz;
    PPzt(:,i) = PPz;
    DOPzt(:,i) = DOPz;
    DOCzt(:,i) = DOCz;
    DICzt(:,i) = DICz;
    O2zt(:,i) = O2z;
    POCzt(:,i) = POCz;
    CO2zt(:,i) = CO2z;
    O2_sat_relt(:,i) = O2_sat_rel;
    O2_sat_abst(:,i) = O2_sat_abs;
    T_sedt(:,i,1) = T_sedim;
    T_sedt(:,i,2) = T_sedim_m;
    T_sedt(:,i,3) = T_sedim_b;
    T_sedt(:,i,4) = T_sedim_t;
    BODzt = 0; %for compatibility with the other code
    %Fokema
    
    DOCzt1(:,i) = DOCz1./DOCz; %Fokema-model DOC subpool 1
    DOCzt2(:,i) = DOCz2./DOCz; %Fokema-model DOC subpool 2
    DOCzt3(:,i) = DOCz3./DOCz; %Fokema-model DOC subpool 3
    DOCtfrac(:,i,1) = DOCz1; %Fokema-model subpool 1 fraction
    DOCtfrac(:,i,2) = DOCz2; %Fokema-model subpool 2 fraction
    DOCtfrac(:,i,3) = DOCz3; %Fokema-model subpool 3 fraction
    Daily_BB1t(:,i) = Daily_BB1; %Fokema-model subpool 1 daily bacterial decomposition
    Daily_BB2t(:,i) = Daily_BB2; %Fokema-model subpool 2 daily bacterial decomposition
    Daily_BB3t(:,i) = Daily_BB3; %Fokema-model subpool 3 daily bacterial decomposition
    Daily_PBt(:,i) = Daily_PB; %Fokema-model daily photobleaching
    
    Qzt_sed(:,i) = Qz_sed./(60*60*24*dt); %(J m-2 day-1) -> (W m-2)
    lambdazt(:,i) = lambdaz_wtot_avg;
    
    surfaceflux(1,i) = surfflux; %Carbon dioxide surface flux
    CO2_eqt(1,i) = CO2_eq;       %Carbon dioxide equilibrium concentration
    K0t(:,i) = K0;               %Dissolved carbon doxide solubility coefficient
    CO2_ppmt(:,i) = CO2_ppm;
    
    O2fluxt(1,i) = O2flux;       %Oxygen surface flux
    O2_eqt(1,i) = O2_eq;         %Oxygen equilibrium concentration
    K0_O2t(1,i) = K0_O2;         %Dissolved oxygen solubility coefficient
    dO2Chlt(:,i) = dO2_Chl;      %DO production by photosynthesis
    dO2SODt(:,i) = dO2_SOD;      %DO consumption by sediment degradation
    dO2DOCt(:,i) = dO2_DOC;      %DO consumption by DOC degradation
    pHt(:,i) = -log10(Hplusz); %pH
    
    POC1_frac = POCz1./POCz; %Fraction of autochthonous POC
    POC2_frac = POCz2./POCz; %Fraction of allochthonous POC
    POC1tfrac(:,i) = POC1_frac;
    
    testi3t(:,i,1) = Growth_bioz;
    testi3t(:,i,2) = Loss_bioz;
    testi3t(:,i,3) = lambdaz_wtot;
    testi3t(:,i,4) = lambdaz_wtot_avg;
    testi3t(:,i,5) = PAR_z;
    testi3t(:,i,6) = POC2_frac;
    
    P3zt_sed(:,i,1) = Pdz_store; %diss. P conc. in sediment pore water (mg m-3)
    P3zt_sed(:,i,2) = Psz_store; %P conc. in inorganic sediment particles (mg kg-1 dry w.)
    P3zt_sed(:,i,3) = Chlsz_store; %Chl conc. in organic sediment particles (mg kg-1 dry w.)
    P3zt_sed(:,i,4) = F_IM; %VOLUME fraction of inorganic particles of total dry sediment
    P3zt_sed(:,i,5) = Sedimentation_counter; %H_netsed_inorg; %Sedimentation (m/day) of inorganic particles of total dry sediment
    P3zt_sed(:,i,6) = Resuspension_counter; %H_netsed_org; %Sedimentation (m/day) of organic particles of total dry sediment
    P3zt_sed(:,i,7) = H_netsed_inorg./H_totsed; %(monitoring variables)
    P3zt_sed(:,i,8) = POC_sedim_counter;
    P3zt_sed(:,i,9) = F_NLOM;
    P3zt_sed(:,i,10) = F_LOM;
    P3zt_sed(:,i,11) = POC_resusp_counter;
    
    P3zt_sed_sc(:,i,1) = dPW_up; %(mg m-3 day-1) change in Pz due to exchange with pore water
    P3zt_sed_sc(:,i,2) = dPP; %(mg m-3 day-1)
    P3zt_sed_sc(:,i,3) = dChl_res; %(mg m-3 day-1)
    
    His(1,i) = Hi;
    His(2,i) = (rho_fw/rho_snow)*WEQs;
    His(3,i) = Hsi;
    His(4,i) = Tice;
    His(5,i) = Wt(i,3);
    His(6,i) = rho_snow;
    His(7,i) = IceIndicator;
    His(8,i) = HFrazil; %NEW!!!
    
    %Original MixStat matrix in v.1.2.1b
    
    %MixStat(1,i) = Iflw_S;
    %MixStat(2,i) = Iflw_TP;
    %MixStat(3,i) = sum(Sz.*Vz);
    %MixStat(4,i) = Growth_bioz(1);%mean(Growth_bioz(1:4)); %Obs! changed to apply to layers 1-4 only
    %MixStat(5,i) = Loss_bioz(1);%mean(Loss_bioz(1:4)); %Obs! changed to apply to layers 1-4 only
    %%MixStat(6,i) = Iflw;
    %MixStat(7:10,i) = NaN;
    
    % MixStat matrix from v.1.2 for figure output purposes
    
    MixStat(1,i) = sum(1e-6*(dPP).*Vz);%Iflw_S;
    MixStat(2,i) = sum(1e-6*(dPW_up).*Vz);%Iflw_TP;
    MixStat(3,i) = sum(1e-6*(dPOC_res/C_per_P_alloc).*Vz);%lambdaz_wtot(2);%Iflw_DOC;
    MixStat(4,i) = mean(Growth_bioz); %Only for chlorophyll group 1 (a)
    MixStat(5,i) = mean(Loss_bioz);  %Only for chlorophyll group 1 (a)
    MixStat(6,i) = Iflw;
    if (IceIndicator == 1)
        MixStat(7:11,i) = NaN;
    else
        dum=interp1(zz,Pz,[0:0.1:4]);
        MixStat(7,i) = mean(dum); %diss-P conc. 0-4m in ice-free period
        
        dum=interp1(zz,Chlz,[0:0.1:4]);
        MixStat(8,i) = mean(dum); %Chla conc. 0-4m in ice-free period
        
        dum=interp1(zz,PPz,[0:0.1:4]);
        MixStat(9,i) = mean(dum); %particulate inorg. P conc. 0-4m in ice-free period
        
        dum=interp1(zz,DOPz,[0:0.1:4]);
        MixStat(10,i) = mean(dum); %dissolved organic P conc. 0-4m in ice-free period
        
        dum=interp1(zz,Sz,[0:0.1:4]);
        MixStat(11,i) = mean(dum); %particulate matter conc. 0-4m in ice-free period
    end
    
    MixStat(12,i) = TCz; %pycnocline depth
    
    MixStat(13,i) = 1e-6*Iflw*Iflw_TP; %total P inflow (kg day-1)
    if (Iflw>Vz(1))
        disp('Large inflow!!')
    end
    MixStat(14,i) = 1e-6*Iflw*(Pz(1)+PPz(1)+DOPz(1)+(Chlz(1)+Cz(1))/Y_cp+POCz1(1)/(C_per_P_pp)+POCz2(1)/C_per_P_alloc); %total P outflow (kg day-1) Qs;
    MixStat(15,i) = sum(1e-6*Vz.*(delPP_inorg + delC_org + delC_org2 + delC_porg./C_P_org_sed)); %total P sink due to sedimentation (kg day-1) Ql;
    MixStat(16,i) = sum(1e-6*(dPP+dPW_up+dPOC_res./C_P_org_sed).*Vz); %Internal P loading (kg day-1, excluding Chla)
    MixStat(17,i) = sum(1e-6*dChl_res.*Vz); %Internal Chla loading (kg day-1)(resuspension 50/50 between the two groups)
    MixStat(18,i)= sum(1e-6*Vz.*((Pz+PPz+DOPz+(Chlz+Cz)/Y_cp+POCz1/(C_per_P_pp)+POCz2/C_per_P_alloc) - TP0)); %Net P change kg
    MixStat(19,i)= sum(1e-6*((dPP+dPW_up-delPP_inorg+dChl_res-delC_org).*Vz - (1-F_sed_sld)*H_sed*(-diff([Az; 0])).*dPW_down)); %Net P flux from sediment kg
    MixStat(20,i) = 1e-6*Iflw*(Iflw_TP-(Iflw_Chl+Iflw_C)./Y_cp-Iflw_DOP-Fstable*Iflw_S); %total algae-available P inflow (kg day-1)
    if (IceIndicator == 1)
        MixStat(21,i) = NaN;
    else
        dum=interp1(zz,Cz,[0:0.1:4]);
        MixStat(21,i) = mean(dum); %Chl group 2 conc. 0-4m in ice-free period
    end
    MixStat(22,i) = lvlD;%happitassa;;%mean(Growth_bioz_2); %For chlorophyll group 2
    MixStat(23,i) = mean(Loss_bioz_2);  %For chlorophyll group 2
    
    %if(i==65);keyboard;end
    
end; %for i = 1:length(tt)

runtime=toc;

disp(['Total model runtime: ' int2str(floor(runtime/60)) ' min ' int2str(round(mod(runtime,60))) ' s']);
disp(['Reduced SS load due to inconsistencies: '  num2str(round(SS_decr)) ' kg']);

% >>>>>> End of the time loop >>>>>>

% Below are the two functions for calculating tridiagonal matrix Fi for solving the
% 1) diffusion equation (tridiag_DIF_v11), and
% 2) advection-diffusion equation (tridiag_HAD_v11) by fully implicit hybrid exponential numerical scheme,
% based on Dhamotharan et al. 1981,
%'Unsteady one-dimensional settling of suspended sediments', Water Resources Research 17(4), 1125-1132
% code checked by TSA, 16.03.2004


%Inputs:
% Kz    diffusion coefficient at layer interfaces (plus surface) N (N,1)
% U     vertical settling velocity (scalar)
% Vz    layer volumes (N,1)
% Az    layer interface areas (N,1)
% dz    grid size
% dt    time step

%Output:
% Fi    tridiagonal matrix for solving new profile Cz

az = (dt/dz) * [0; Kz] .* (Az ./ Vz);
bz = (dt/dz) * [Kz; 0] .* ([Az(2:end); 0] ./ Vz);
Gi = [-bz (1 + az + bz) -az];

%=== DIFFUSIVE EQUATION ===
function[Fi]=tridiag_DIF_v11(Kz,Vz,Az,dz,dt)

Nz=length(Vz); %number of grid points/layers

% Linearized heat conservation equation matrix (diffusion only)
az = (dt/dz) * Kz .* (Az ./ Vz);                                        %coefficient for i-1
cz = (dt/dz) * [Kz(2:end); NaN] .* ([Az(2:end); NaN] ./ Vz);            %coefficient for i+1
bz = 1 + az + cz;                                                       %coefficient for i+1
%Boundary conditions, surface

az(1) = 0;
%cz(1) remains unchanged
bz(1)= 1 + az(1) + cz(1);


%Boundary conditions, bottom

%az(end) remains unchanged
cz(end) = 0;
bz(end) = 1 + az(end) + cz(end);

Gi = [-cz bz -az];
Fi = spdiags(Gi,-1:1,Nz,Nz)';
%end of function


%=== ADVECTIVE-DIFFUSIVE EQUATION ===
function[Fi]=tridiag_HAD_v11(Kz,U,Vz,Az,dz,dt)

if (U<0)
    error('only positive (downward) velocities allowed')
end

if (U==0)
    U=eps; %set Vz next to nothing (=2.2204e-016) in order to avoid division by zero
end

Nz=length(Vz); %number of grid points/layers

theta=U*(dt/dz);

az = theta.*(1 + (1./(exp( (U*Vz)./(Kz.*Az) ) - 1)));                   %coefficient for i-1
cz = theta./(exp( (U*Vz)./([Kz(2:end); NaN].*[Az(2:end); NaN]) ) - 1);  %coefficient for i+1
bz = 1 + az + cz;                                                       %coefficient for i

%Boundary conditions, surface

az(1) = 0;
%cz(1) remains unchanged
bz(1) = 1 + theta + cz(1);

%Boundary conditions, bottom

%az(end) remains unchanged
cz(end) = 0;
bz(end) = 1 + az(end);

Gi = [-cz bz -az];
Fi = spdiags(Gi,-1:1,Nz,Nz)';
%end of function


function [Pdiss, Pfpart]=Ppart(vf,TIP,Psat,Fmax,rho_sed,Fstable)
% Function for calculating the partitioning between
% dissolved and inorganic particle bound phosphorus.
% Based on Langmuir isotherm approach
%vf:    volume fraction of suspended inorganic matter (m3 m-3); S/rho_sed OR (1-porosity)
%TIP:   Total inorganic phosphorus (mg m-3)
%Psat, mg m-3 - Langmuir half-saturation parameter
%Fmax, mg kg-1  - Langmuir scaling parameter
%rho_sed, kg m-3 - Density of dry inorganic sediment mass
%Fstable, mg kg-1 - Inactive P conc. in inorg. particles

N=length(TIP);
Pdiss=NaN*ones(N,1);

for w=1:N
    a = vf(w)-1;
    b = TIP(w) + (vf(w)-1)*Psat - vf(w)*rho_sed*(Fmax+Fstable);
    c = Psat*TIP(w) - vf(w)*rho_sed*Fstable*Psat ;
    try
        Pdiss(w) = max(real(roots([a b c])));
    catch
        keyboard
    end
end


%NEW!!!! Treshold value for numerical stability (added 020707):
%truncate negative values
cutinx=find(Pdiss < 0);
if (isempty(cutinx)==0)
    Pdiss(cutinx)=0;
    disp('NOTE: Pdiss < 0, values truncated')
end

%truncate too high values
cutinx=find(Pdiss > (TIP - Fstable*rho_sed*vf));
if (isempty(cutinx)==0)
    Pdiss(cutinx)=(TIP(cutinx) - Fstable*rho_sed*vf(cutinx));
    disp('NOTE: Pdiss > (TIP - Fstable*rho_sed*vf), values truncated')
end

Pfpart = (TIP - (1-vf).*Pdiss)./(rho_sed*vf); %inorg. P conc. in sediment particles(mg kg-1 dry w.)
%end of function

