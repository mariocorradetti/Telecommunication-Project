%Rain height model for prediction methods ITU-R P.839 

LON_hR = dlmread('rainHeight/Lon.txt');
lon_hR = LON_hR(1,:)'; % from 0 to 360
LAT_hR = dlmread('rainHeight/Lat.txt');
lat_hR = LAT_hR(:,1);  % from 90 to -90

% Mean annual 0-celsius isotherm height [km] 
h0_map  = dlmread('rainHeight/h0.txt'); % (LAT, LON)
% Shift map s.t. Longitude goes from -180 to 180
[lon_hR, h0_map] = MapData_centerGreenwich(lon_hR, h0_map);
[LON_hR, LAT_hR] = meshgrid(lon_hR, lat_hR);

% Mean annual rain height [km]
hR_map = h0_map + 0.36;

Map.LON_hR = LON_hR;
Map.lon_hR = lon_hR;
Map.LAT_hR = LAT_hR;
Map.lat_hR = lat_hR;
Map.hR_map = hR_map;

% Characterisitcs of precipitation for prediction methods ITU-R P.837 

LON_R001 = dlmread('rainfallRate/LON_R001.TXT');
lon_R001 = LON_R001(1,:)'; % from -180 to 180
LAT_R001 = dlmread('rainfallRate/LAT_R001.TXT');
lat_R001 = LAT_R001(:,1); % from 90 to -90

%Rainfall rate exceeded for 0.01% of an average year [mm/h] 
R001_map  = dlmread('rainfallRate/R001.TXT'); % (LAT, LON)

Map.LON_R001 = LON_R001;
Map.lon_R001 = lon_R001;
Map.LAT_R001 = LAT_R001;
Map.lat_R001 = lat_R001;
Map.R001_map = R001_map;

% Attenuations due to cloudsand fog: ITU-R P.840 %%%%%%%%%%%

LON_Lred = dlmread('cloudWaterContent/ESALON1dot125.txt');
lon_Lred = LON_Lred(1,:)'; % from 0 to 360
LAT_Lred = dlmread('cloudWaterContent/ESALAT1dot125.txt');
lat_Lred = LAT_Lred(:,1);  % from 90 to -90

% Reduced columnar cloud liquid water content exceeded for 0.1% of an
% average year
Lred01_map  = dlmread('cloudWaterContent/Lred_1_v4.txt'); % (LAT, LON)

[lon_Lred, Lred01_map] = MapData_centerGreenwich(lon_Lred, Lred01_map);
[LON_Lred, LAT_Lred] = meshgrid(lon_Lred, lat_Lred);

Map.LON_Lred   = LON_Lred;
Map.lon_Lred   = lon_Lred;
Map.LAT_Lred   = LAT_Lred;
Map.lat_Lred   = lat_Lred;
Map.Lred01_map = Lred01_map;

% Attenuation by atmospheric gases and related effects ITU-R P.676: Spectroscopic data for oxygen and water vapour attenuation

Table.ox = dlmread('atm_properties/oxygen_spectrData.txt');
Table.wv = dlmread('atm_properties/waterVapor_spectrData.txt');
Table.wv_height = dlmread('atm_properties/waterVapor_heightData.txt');

% Water vapour ITU-R P.836 : water vapour density at earth surface from the maps

Map.LON_rho = dlmread('atm_properties/LON1dot125.txt');
Map.lon_rho = Map.LON_rho(1,:)'; % from 0 to 360
Map.LAT_rho = dlmread('atm_properties/LAT1dot125.txt');
Map.lat_rho = Map.LAT_rho(:,1);  % from 90 to -90
Map.rho1_map  = dlmread('atm_properties/RHO_1_v4.txt'); % (LAT, LON)
[Map.lon_rho, Map.rho1_map] = MapData_centerGreenwich(Map.lon_rho, Map.rho1_map);
[Map.LON_rho,Map.LAT_rho] = meshgrid(Map.lon_rho,Map.lat_rho);

% Refractive index ITU-R P.453 : wet term of the radio refractivity from the maps

Map.LON_Nwet = dlmread('atm_properties/LON_N.txt');
Map.lon_Nwet = Map.LON_Nwet(1,:)'; % from -180 to 180
Map.LAT_Nwet = dlmread('atm_properties/LAT_N.txt');
Map.lat_Nwet = Map.LAT_Nwet(:,1);  % from 90 to -90
Map.Nwet1_map  = dlmread('atm_properties/NWET_Annual_1.txt'); % (LAT, LON)
