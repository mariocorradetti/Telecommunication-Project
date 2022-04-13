%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                 %
%                 TELECOMMUNICATION SYSTEMS PROJECT               %
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

% LOAD ATTENUATION DATA
load ../savedData/attenuationMap.mat
A_gas = Agas;
A_rain = Arain001;
A_cloud = Acloud001;

run('System')
run('GroundStation.m')


%% ATTENUATION MAP WITH 5 deg ELEVATION
% ITU-R P618-13 : overall attenuation for a give probability (1%)
A_tot = Agas + sqrt((Arain001+Acloud001).^2 ); %(60)

% Plot 
h = figure;
plotMap(lon_vec,lat_vec, A_tot)
hcb = colorbar;
set(get(hcb,'Title'),'String','$A_\mathrm{1\%}$ [dB]', 'Interpreter', 'latex', 'FontSize', fontSize);
title('Overall atmospheric attenuation exceeded for 1\% of an average year ')
set(h, 'Units', 'normalized', 'OuterPosition', [.1 .3 .7 .6]);

%% ATTENUATION AT SELECTED G/S SITE WITH DIFFERENT ELEVATION ANGLES                  

th_vec = deg2rad([rad2deg(th_el):90]); % vector of elevation angle
% Plot
h = figure;
set(h, 'Units', 'normalized', 'OuterPosition', [.1 .3 .7 .6]);
grid on
grid minor
hold on

for j = 1:length(th_vec)
    % Selecting elevation angle
    th = th_vec(j);
    % Attenuation evaluation 
    A_cloud = cloudAtt001(f, th, GS1, Map);
    A_rain = rainAtt001(f, th, GS1, tau, Map);
    [gam_o, h_o] = gas_specAtt(f,'oxygen',GS1,Map,Table);
    [gam_w, h_w] = gas_specAtt(f,'waterVapour',GS1,Map,Table);
    A_gas = inclined_path(gam_o,gam_w,h_o,h_w,th,GS1.h,h2,'earth-space');
    % Overall attenuation 
    A_tot_GS(j) = A_gas + sqrt((A_rain+A_cloud).^2 );

end

plot(rad2deg(th_vec), A_tot_GS(:),'LineWidth', 1)    

xlabel('Elevation angle [deg]')
ylabel('Attenuation [dB]')
title('Atmospheric attenuation at New Norcia ground station')

%% SAVE DATA
save('../savedData/attenuationMap.mat', 'A_tot', '-append')
save('../savedData/attenuationMap.mat', 'A_tot_GS', '-append')
save('../savedData/attenuationMap.mat', 'th_vec', '-append')


