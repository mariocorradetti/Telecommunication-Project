
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                 %
%               TELECOMMUNICATION SYSTEMS PROJECT                 %
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
fontSize = 15;
set(0,'DefaultAxesFontSize',fontSize);

addpath('..\functions')
addpath('..\functions\time\')



% Set useful data
data.mu = astroConstants(13);  % Earth’s gravitational parameter [km^3/s^2]
data.we = deg2rad(15.04/3600); % Earth’s rotation velocity [km/s]
data.J2 = astroConstants(9);   % Gravitatonal field constant of the Earth  
data.Re = astroConstants(23);  % Earth's mean radius [km]
data.th_g0 = 0;                % Longitude of Greenwich meridian at initial time [rad]

% Given orbit 
a = 0.7150e4;   % Semi-major axis [km]
e = 0*0.0054;     % Eccentricity [-]
i = deg2rad(47.9875); % Inclination [rad]
OM = deg2rad(260);  % Our choice [rad]
om = deg2rad(-180); % Our choice [rad]
th = deg2rad(50);   % Our choice [rad]

% -------------------------------------------------------------------------
% Set time vectors for ground tracks
T_orb = 2*pi*sqrt(a^3/data.mu); % s/c orbital period [s]
T_day = 24*3600;                % One day [s]


%% NOMINAL GROUNDTRACK
t0 = 0;     % Initial time [s]
norb = 1;   % Number of s/c revolutions for ground track

tspan1 = linspace(t0, norb*T_orb, 10000);   % time span over one orbit
tspan2 = linspace(t0, T_day, 10000);        % time span over one day

% Compute initial state vector
kep = [a,e,i,OM,om,th];             % keplerian elements
[rr0,vv0] = kep2car(kep, data.mu);  % compute initial position and velocity vectors

% Ground track with J2 effect
code = 1 ; % takes into account Earth oblateness (J2 effect)
[lonp1,latp1] = groundtrack(rr0, vv0, tspan1, data, code); % gt over 1 orbit
[lonp2,latp2] = groundtrack(rr0, vv0, tspan2, data, code); % gt over 1 day

figure()
subplot(2,1,1)
print_GT(lonp1, latp1, 'Groundtrack of the orbit')
subplot(2,1,2)
print_GT(lonp2, latp2, '1 day with J2')


%% REPEATING GROUNDTRACK 
k = 14; % Revolutions of the satellite 
m = 1;  % Rotations of the planet 

% Compute new semi-major axis considering oblateness effects 
const = -3/2*sqrt(data.mu)* data.J2* data.Re^2/(1-e^2)^2;
fun = @(a) m/k-(data.we-const/a^(7/2)*cos(i))/(sqrt(data.mu/a^3)+const/a^(7/2)*(5/2*sin(i)^2-2)+const/a^(7/2)*(1-3/2*(sin(i))^2));
a_rep = fzero(fun, a);

% new orbital elements
kep = [a_rep, e, i, OM, om, th];
[rr0,vv0] = kep2car(kep, data.mu);

% New orbital period 
T_orb = 2*pi*sqrt(a_rep^3/ data.mu);      
tspan1 = linspace(t0, norb*T_orb, 1000);  

% Ground track with J2
code = 1;
[lonp1,latp1] = groundtrack(rr0, vv0, tspan1, data, code);
[lonp2,latp2] = groundtrack(rr0, vv0, tspan2, data, code);

figure()
print_GT(lonp2, latp2, 'Repeating groundtrack: 14 satellite revolutions')
hold on
% plot ground station locations
scatter(-3-57/60-5.7/3600, 40+26/60+33.23/3600,'filled', 'MarkerFaceColor','r','MarkerEdgeColor','w');
scatter(-25-8/60-8.6/3600, 36+59/60+50.1/360,'filled', 'MarkerFaceColor','r','MarkerEdgeColor','w');
scatter(115-53/60-6.58/3600, -31-48/60-9.08/3600,'filled', 'MarkerFaceColor','r','MarkerEdgeColor','w');



%% VISIBILITY WINDOW
% Orbit propagation 
% data.y0 = [rr0; vv0];             
% code = 1;       % consider Earth oblateness
% [t, y] = ode113(@(t,y) ode_2body_per(t,y,data,code), tspan1, data.y0);
% r = y(:,1:3);

% Ground station coordinates ----------------------------------------------
% Villafranca - Spain
lat1 = (40+26/60+33.23/3600);   % [deg]
lon1 = (-3-57/60-5.7/3600);     % [deg]
h1 = 664.80/1000;               % [km] 

% Perth
lat2 = (-31-2/60-54/3600);      % [deg]
lon2 = (116+11/60+28/3600);     % [deg]
h2 = 22.16/1000;                % [km]
% -------------------------------------------------------------------------

lat_gs_vect = [lat1, lat2];
lon_gs_vect = [lon1, lon2];
h_gs_vect = [h1,h2];

% Meshgrid for elevation angle plot
lon_vec = -180:0.25:180;
lat_vec = -90:0.25:90;
[LON_VEC, LAT_VEC] = meshgrid(lon_vec, lat_vec);


theta_map = zeros(length(lat_vec), length(lon_vec));
gs_map = zeros(size(theta_map)); 

for k = 1:length(lon_gs_vect)
    lat_gs = lat_gs_vect(k);
    lon_gs = lon_gs_vect(k);
    
    % Distance from Earth's center of the GS [km]
    R_local = data.Re + h_gs_vect(k);

    r_GS = R_local * [cosd(lat_gs)*cosd(lon_gs), cosd(lat_gs)*sind(lon_gs), sind(lat_gs)]';

    for m = 1:length(lat_vec)
        lat_sc = lat_vec(m);

        for n = 1:length(lon_vec)
            lon_sc = lon_vec(n);

            r_SC = a_rep * [cosd(lat_sc)*cosd(lon_sc), cosd(lat_sc)*sind(lon_sc), sind(lat_sc)]';
            delta_r = r_SC - r_GS;

            compl_theta = acosd( dot(delta_r, r_GS) / (norm(delta_r)*norm(r_GS)) ); % Complementary angle of elevation angle [deg]            
            if compl_theta <= 90 %theta_map(m,n)
                theta_map(m,n) = 90-compl_theta;
                gs_map(m,n) = k;
            end

        end
    end

% Assuming a 5 deg elevation limit
elevation_limit = 5; % [deg]
theta_map(theta_map < elevation_limit) = NaN;
gs_map(isnan(theta_map)) = 0; % Set to zero instead of NaN for computing the visibility time
alphaData = gs_map>0;


h = figure;
set(h, 'Units', 'normalized');
plotMap(lon_vec,lat_vec, theta_map, alphaData)
grid off
colormap(gca, jet)
hcb = colorbar;
caxis([min(min(theta_map)) max(max(theta_map))])
set(get(hcb,'Title'),'String','$\theta$ [deg]', 'Interpreter', 'latex');
title('Elevation angle')
hold on
plot(lonp2, latp2, 'r', 'LineWidth', 1, 'HandleVisibility', 'off');
end



%% VISIBILITY TIME
clear theta_visibility_time

% Select the ground station -----------------------------------------------
% % Villafranca - Spain
% lat_gs = (40+26/60+33.23/3600);
% lon_gs = (-3-57/60-5.7/3600);
% h_gs = 664.80; % m

% % Santa Maria - Portogallo
% lat_gs = (36+59/60+50.1/3600);
% lon_gs = (-25-8/60-8.6/3600);
% h_gs = 276;
% % 
% New Norcia
lat_gs = (-31-2/60-54/3600);
lon_gs = (116+11/60+28/3600);
h_gs = 22.16; % m
% -------------------------------------------------------------------------

% Distance from Earth's center of the GS [km]
R_local = data.Re + h_gs;

r_GS = R_local * [cosd(lat_gs)*cosd(lon_gs), cosd(lat_gs)*sind(lon_gs), sind(lat_gs)]';

for i = 1:length(lonp2)
    lat_sc = latp2(i);
    lon_sc = lonp2(i);

    r_SC = a_rep * [cosd(lat_sc)*cosd(lon_sc), cosd(lat_sc)*sind(lon_sc), sind(lat_sc)]';
    delta_r = r_SC - r_GS;

    compl_theta = acosd( dot(delta_r, r_GS) / (norm(delta_r)*norm(r_GS)) ); % Complementary angle of elevation angle [deg]            
    if compl_theta <= 90 %theta_map(m,n)
        theta_visibility_time(i) = 90-compl_theta;
    end
end

% Assuming a 5 deg elevation limit
elevation_limit = 5; % [deg]
theta_visibility_time(theta_visibility_time < elevation_limit) = 0;
time_counter = 0;

for i = 1:length(theta_visibility_time)
    if theta_visibility_time(i) == 0
        lon_plot(i) = NaN;
        lat_plot(i) = NaN;
        time_vis_vect(i,1) = 0;
    else
        lon_plot(i) = lonp2(i);
        lat_plot(i) = latp2(i);
        time_counter = time_counter +1;
        time_vis_vect(i,1) = 1;        
    end
end

dt = tspan2(2)-tspan2(1);
time_vis = time_counter*dt;
fprintf('\nThe visibility time window above the selected ground station is %f minutes',time_vis/60)

figure()
print_GT(lonp2, latp2, '')
hold on
scatter(lon_gs, lat_gs,'MarkerEdgeColor','k','MarkerFaceColor','w')
plot(lon_plot, lat_plot,'r','LineWidth', 1.5)
title('Visibility region over ground station')

        
%% Data budget
data_rate = 1.06e-3; %[Gbit/s]
data_volume = data_rate*time_vis;
fprintf('\n Total data volume available for download: %f Gbit', data_volume)