function [A1perc] = cloudAtt001(f, th, GS, DigitalMap)
% ITU-R P840-8

LON_Lred   = DigitalMap.LON_Lred  ;
lon_Lred   = DigitalMap.lon_Lred  ;
LAT_Lred   = DigitalMap.LAT_Lred  ;
lat_Lred   = DigitalMap.lat_Lred  ;
Lred01_map = DigitalMap.Lred01_map;


T_red = 273.15; % Equivalent cloud liquid water temperature [K]

theta = 300 / T_red; % Normalized temperature
eps0 = 77.66 + 103.3*(theta-1);
eps1 = 0.0671 * eps0;
eps2 = 3.52;
f_p = 20.20 - 146*(theta-1) + 316*(theta-1)^2; % Principal relaxation frequency [GHz]
f_s = 39.8*f_p; % Secondary relaxation frequency [GHz]
% Complex dielectric permittivity of water:
eps_prime = (eps0-eps1) / (1+(f/f_p)^2) + (eps1-eps2) / (1+(f/f_s)^2) + eps2;
eps_second = f*(eps0-eps1) / (f_p*(1+(f/f_p)^2)) + f*(eps1-eps2) / (f_s*(1+(f/f_s)^2));

eta = (2+eps_prime) / eps_second;

% Cloud liquid water specific attenuation coefficient
K_l = 0.819*f / (eps_second*(1+eta^2)); % [(dB/km)/(g/m^3)] (2)

% Reduced columnar cloud liquid water content exceeded for 0.1% of an
% average year
Lred_1perc = interp2(LON_Lred, LAT_Lred, Lred01_map, GS.lon, GS.lat); % [kg/m^2]

% Cloud attenuation exceeded for 0.1% of an average year
A1perc = Lred_1perc * K_l / sin(th);


end

