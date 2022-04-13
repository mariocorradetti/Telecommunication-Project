function [A_p] = rainAtt001(f, th, GS, tau, DigitalMap)
% f :  frequency [GHz]
% th : path elevation angles [deg]
% tau : polarization tilt angle relative to the horizontal
 
LON_hR = DigitalMap.LON_hR;
lon_hR = DigitalMap.lon_hR;
LAT_hR = DigitalMap.LAT_hR;
lat_hR = DigitalMap.lat_hR;
hR_map = DigitalMap.hR_map;

LON_R001 = DigitalMap.LON_R001;
lon_R001 = DigitalMap.lon_R001;
LAT_R001 = DigitalMap.LAT_R001;
lat_R001 = DigitalMap.lat_R001;
R001_map = DigitalMap.R001_map;

% Procedure used in ITU-R P618-13

% STEP 1: estimate RAIN HEIGHT, using ITU-R P839
hR = interp2(LON_hR, LAT_hR, hR_map, GS.lon, GS.lat);
if (hR-GS.h) <= 0
    % i.e. zero rain attenuation
    A_p = 0;
    return
end

% STEP 2: compute slant path length: L_s [km] (below the rain height)
if th >= deg2rad(5)
    L_s = (hR-GS.h) / sin(th);
else
    R_eff = 8500; % Effective Earth's radius [km]
    L_s = ( 2*(hR-GS.h) ) / ( sqrt((sin(th))^2 + (2*(hR-GS.h))/R_eff) + sin(th));    
end

% STEP 3: compute horizontal projection of L_s
L_G = L_s * cos(th); % [km]

% STEP 4: estimate rainfall rate, using ITU-R P.837
% R001 rainfall rate exceeded for 0.01% of an average year with an
% integration time of 1 min
R001 = interp2(LON_R001, LAT_R001, R001_map, GS.lon, GS.lat);
if R001 == 0
    % i.e. zero rain attenuation
    A_p = 0;
    return
end

% STEP 5: compute specific rain attenuation : gammaR [dB/km]
% coefficient take form ITU-R P838
[kH, kV, alphaH, alphaV] = attenuationCoeffs(f);

k = (kH + kV + (kH-kV)*(cos(th))^2*cos(2*tau)) / 2;
alpha = (kH*alphaH + kV*alphaV + (kH*alphaH-kV*alphaV)*(cos(th))^2*cos(2*tau)) / (2*k);

gammaR = k * R001^alpha; % [dB/km]

% STEP 6: compute horizobtal reduction factor for 0.01% of the time : r001
r001 = ( 1 + 0.78*sqrt((L_G*gammaR)/f) -0.38*(1-exp(-2*L_G)) )^-1;

% STEP 7: compute vertical adjustment factor for 0.01% of the time: nu001
zeta = atand( (hR-GS.h)/(L_G*r001) ); % [deg]

if zeta > rad2deg(th)
    L_R = L_G*r001 / cos(th); % [km]
else
    L_R = (hR-GS.h) / sin(th); % [km]
end

if abs(GS.lat) < 36
    chi = 36 - abs(GS.lat); % [deg]
else
    chi = 0; % [deg]
end

nu001 = ( 1 + sqrt(sin(th)) * (31* (1-exp(-(th/(1+chi)))) * sqrt(L_R*gammaR)/(f^2) - 0.45))^-1;

% STEP 8: compute effective path length
L_E = L_R * nu001; % [km]

% STEP 9: compute attenuation exceeded for 0.01% of an average year [dB]
A001 = gammaR * L_E; % [dB]

% STEP 10: extrapolation for other excedance percentages
p = 1; % extrapolate A_1%

if p >= 1 || abs(GS.lat) >= 36
    beta = 0;
elseif p < 1 && abs(GS.lat) < 36 && th >= deg2rad(25)
    beta = -0.005*(abs(GS.lat) - 36);
else
    beta = -0.005*(abs(GS.lat) - 36) + 1.8 - 4.25*sin(th);
end

A_p = A001 * (p/0.01) ^ -( .655 + .033*log(p) -0.045*log(A001) - beta*(1-p)*sin(th));



end

