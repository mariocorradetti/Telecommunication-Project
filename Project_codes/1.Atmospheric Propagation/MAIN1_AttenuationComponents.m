%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                 %
%              TELECOMMUNICATION SYSTEMS PROJECT                  %
%                    Academic year 2020/2021                      %
%                                                                 %
%                     Balossi Claudia                             %   
%                     Corradetti Mario                            %
%                     Donato Giuseppe                             %
%                                                                 %
%                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc
close all
clear all 

set(0,'defaultTextInterpreter','latex')
fontSize = 20;
set(0,'DefaultAxesFontSize', fontSize);

addpath('../functions')
addpath('../savedData')

load topo.mat % Load Earth's topography [m] implemented in Matlab with toolbox

run('System')


% Propagation data and prediction methods required for the design 
% of Earth-space telecommunication systems      
% Recommendation : ITU-R P.618-13 (12/2017)

%% RAIN ATTENUATION
% Rain height model for prediction methods: ITU-R P.839 (P.839-4)

LON_hR = dlmread('rainHeight/Lon.txt');
lon_hR = LON_hR(1,:)'; % from 0 to 360
LAT_hR = dlmread('rainHeight/Lat.txt');
lat_hR = LAT_hR(:,1);  % from 90 to -90
% Mean annual 0°C isotherm height [km] 
h0_map  = dlmread('rainHeight/h0.txt'); %(LAT,LON)
% Shift map so the Longitude goes from -180 to 180
[lon_hR, h0_map] = MapData_centerGreenwich(lon_hR, h0_map);
[LON_hR, LAT_hR] = meshgrid(lon_hR, lat_hR);
% Mean annual rain height [km]
hR_map = h0_map + 0.36;

% Plot 
h = figure;
plotMap(lon_hR,lat_hR, hR_map)
hcb = colorbar;
set(get(hcb,'Title'),'String','$h_R$ [km]', 'Interpreter', 'latex', 'FontSize', fontSize);
title('Mean annual rain height')
set(h, 'Units', 'normalized', 'OuterPosition', [.1 .3 .7 .6]);


%Save data in a structure
Map.LON_hR = LON_hR;
Map.lon_hR = lon_hR;
Map.LAT_hR = LAT_hR;
Map.lat_hR = lat_hR;
Map.hR_map = hR_map;

% Characteristics of precipitation for propagation modelling : ITU-R P.837
% (P.837-7)

LON_R001 = dlmread('rainfallRate/LON_R001.TXT');
lon_R001 = LON_R001(1,:)'; % from -180 to 180
LAT_R001 = dlmread('rainfallRate/LAT_R001.TXT');
lat_R001 = LAT_R001(:,1); % from 90 to -90
%Rainfall rate exceeded for 0.01% of an average year [mm/h] 
R001  = dlmread('rainfallRate/R001.TXT'); % (LAT, LON)

% Plot
h = figure;
plotMap(lon_R001,lat_R001, R001)
hcb = colorbar;
set(get(hcb,'Title'),'String','$R_\mathrm{0.01\%}$ [mm/h]', 'Interpreter', 'latex', 'FontSize', fontSize);
title('Rainfall rate exceeded for 0.01\% of an average year')
set(h, 'Units', 'normalized', 'OuterPosition', [.1 .3 .7 .6]);

% Save data
Map.LON_R001 = LON_R001;
Map.lon_R001 = lon_R001;
Map.LAT_R001 = LAT_R001;
Map.lat_R001 = lat_R001;
Map.R001_map = R001;

%% GASEOUS ATTENUATION
% Attenuation by atmospheric gases and related effects ITU-R P.676
% (P.676-12)

Table.ox = dlmread('atm_properties/oxygen_spectrData.txt');
Table.wv = dlmread('atm_properties/waterVapor_spectrData.txt');
Table.wv_height = dlmread('atm_properties/waterVapor_heightData.txt');

% Water vapour: surface density and total columnar content ITU-R P.836
% (P.836-6)
% IMPORT DATA: water vapour density at Earth surface  
Map.LON_rho = dlmread('atm_properties/LON1dot125.txt');
Map.lon_rho = Map.LON_rho(1,:)'; % from 0 to 360
Map.LAT_rho = dlmread('atm_properties/LAT1dot125.txt');
Map.lat_rho = Map.LAT_rho(:,1);  % from 90 to -90
Map.rho1_map  = dlmread('atm_properties/RHO_1_v4.txt'); % (LAT, LON)
[Map.lon_rho, Map.rho1_map] = MapData_centerGreenwich(Map.lon_rho, Map.rho1_map);
[Map.LON_rho, Map.LAT_rho] = meshgrid(Map.lon_rho, Map.lat_rho);

% Reference Standard armospheres ITU-R P.835 (P.835-6)

% Plot surface temperature and pressure
T_map=zeros(length(Map.lat_rho),length(Map.lon_rho));
p_map=zeros(length(Map.lat_rho),length(Map.lon_rho));
for i = 1:length(Map.lat_rho)
    lat = Map.lat_rho(i);
    for j = 1:length(Map.lon_rho)
        lon = Map.lon_rho(j);
        [T,press]=read_mprof(lat,lon);
        T_map(i,j) = T;
        p_map(i,j) = press/1000;
    end
end
% Plot
h=figure;
plotMap(Map.lon_rho, Map.lat_rho, T_map)
hcb = colorbar;
set(get(hcb,'Title'),'String','T [K]', 'Interpreter', 'latex', 'FontSize', fontSize);
title('Temperature of an average year')
set(h, 'Units', 'normalized', 'OuterPosition', [.1 .3 .7 .6]);

% Plot
h=figure;
plotMap(Map.lon_rho,Map.lat_rho, p_map)
hcb = colorbar;
set(get(hcb,'Title'),'String','p [bar]', 'Interpreter', 'latex', 'FontSize', fontSize);
title('Pressure of an average year')
set(h, 'Units', 'normalized', 'OuterPosition', [.1 .3 .7 .6]);


%% CLOUD ATTENUATION
% Attenuation due to cloud and fog ITU-R P.840 (P.840-8)

LON_Lred = dlmread('cloudWaterContent/ESALON1dot125.txt');
lon_Lred = LON_Lred(1,:)'; % from 0 to 360
LAT_Lred = dlmread('cloudWaterContent/ESALAT1dot125.txt');
lat_Lred = LAT_Lred(:,1);  % from 90 to -90

% Reduced columnar cloud liquid water content exceeded for 0.1% of an
% average year
Lred01_map  = dlmread('cloudWaterContent/Lred_1_v4.txt'); % (LAT, LON)
% Shift map s.t. Longitude goes from -180 to 180
[lon_Lred, Lred01_map] = MapData_centerGreenwich(lon_Lred, Lred01_map);
[LON_Lred, LAT_Lred] = meshgrid(lon_Lred, lat_Lred);

h = figure;
plotMap(lon_Lred,lat_Lred, Lred01_map)
hcb = colorbar;
set(get(hcb,'Title'),'String','$L_\mathrm{1\%}^\mathrm{red}$ [kg/m$^2$]', 'Interpreter', 'latex', 'FontSize', fontSize);
title('Reduced cloud liquid water content exceeded for 1\% of an average year')
set(h, 'Units', 'normalized', 'OuterPosition', [.1 .3 .7 .6]);

Map.LON_Lred   = LON_Lred;
Map.lon_Lred   = lon_Lred;
Map.LAT_Lred   = LAT_Lred;
Map.lat_Lred   = lat_Lred;
Map.Lred01_map = Lred01_map;

%% ATTENUATION EVALUATION
lon_vec = [-180:1.5:180];
lat_vec = [-90:1.5:90];

Arain001  = zeros(length(lat_vec), length(lon_vec));
Acloud001 = zeros(length(lat_vec), length(lon_vec));
Agas       = zeros(length(lat_vec), length(lon_vec));

% Waitbar
wBar = waitbar(0, 'Computing attenuation maps...');
titleHandle = get(findobj(wBar,'Type','axes'),'Title');
set(titleHandle, 'FontSize', 16)

for i = 1:length(lat_vec)
    GS.lat = lat_vec(i);
    waitbar(i/length(lat_vec), wBar, 'Computing attenuation maps...')
    
    for j = 1:length(lon_vec)
        GS.lon = lon_vec(j);
        GS.h = 1e-3 * max(0, topo(max(1,round(90+GS.lat)), max(1,round(180+GS.lon))) ); % [km]
    
        % Rain
        Arain001(i,j) = rainAtt001(f, th_el, GS, tau, Map);
        
        % Clouds
        Acloud001(i,j) = cloudAtt001(f, th_el, GS, Map);
        % Gas
        [gam_o, h_o] = gas_specAtt(f,'oxygen',GS,Map,Table);
        [gam_w, h_w] = gas_specAtt(f,'waterVapour',GS,Map,Table);
        Agas (i,j) = inclined_path(gam_o,gam_w,h_o,h_w,th_el,GS.h,h2,'earth-space');
                       
    end
end
close(wBar)

Agas(end-1,:)=Agas(end-2,:);
Acloud001(end-1,:)=Acloud001(end-2,:);

%% PLOT ATTENUATION
% RAIN
h = figure;
plotMap(lon_vec,lat_vec, Arain001)
hcb = colorbar;
set(get(hcb,'Title'),'String','$A_\mathrm{1\%}$ [dB]', 'Interpreter', 'latex', 'FontSize', fontSize);
title('Rain attenuation exceeded for 1\% of an average year')
set(h, 'Units', 'normalized', 'OuterPosition', [.1 .3 .7 .6]);

% CLOUDS
h = figure;
plotMap(lon_vec,lat_vec, Acloud001)
hcb = colorbar;
set(get(hcb,'Title'),'String','$A_\mathrm{1\%}$ [dB]', 'Interpreter', 'latex', 'FontSize', fontSize);
title('Cloud attenuation exceeded for 1\% of an average year')
set(h, 'Units', 'normalized', 'OuterPosition', [.1 .3 .7 .6]);

% GAS
h = figure;
plotMap(lon_vec,lat_vec, Agas)
hcb = colorbar;
set(get(hcb,'Title'),'String','$A_\mathrm{1\%}$ [dB]', 'Interpreter', 'latex', 'FontSize', fontSize);
title('Gaseous attenuation exceeded for 1\% of an average year')
set(h, 'Units', 'normalized', 'OuterPosition', [.1 .3 .7 .6]);

%% SAVE DATA
save('../savedData/attenuationMap.mat','Map')
save('../savedData/attenuationMap.mat','Table', '-append')
save('../savedData/attenuationMap.mat', 'Arain001', '-append')
save('../savedData/attenuationMap.mat', 'Acloud001', '-append')
save('../savedData/attenuationMap.mat', 'Agas', '-append')
save('../savedData/attenuationMap.mat', 'lon_vec', '-append')
save('../savedData/attenuationMap.mat', 'lat_vec', '-append')
