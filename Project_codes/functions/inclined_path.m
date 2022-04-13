function [A] = inclined_path(gam_o,gam_w,h_o,h_w,phi,h1,h2, type)

if nargin == 5
    type = 'earth-space';
end
if phi >= deg2rad(5)
    % Elevation angles between 5 and 90deg
    
    if strcmp(type,'earth-space')
      % Earth-space paths
      Ao = gam_o*h_o;
      Aw = gam_w*h_w;
      A = (Ao + Aw)/sin(phi); % [dB] % path attenuation is obtained using the cosecant law
    end
    
    if strcmp(type,'inclined')
      % Inclined paths
      h_o_ = h_o * (exp(-h1/h_o) - exp(-h2/h_o)); % [km]
      h_w_ = h_w * (exp(-h1/h_w) - exp(-h2/h_w)); % [km]
      A = gam_o*h_o_ + gam_w*h_w_; % Total zenith attenuation
    end

else
    % Elevation angles between 0 and 5deg
    
    if strcmp(type,'earth-space')
      % Earth-space paths
      disp('see Annex 1')
    end
      
       if strcmp(type,'inclined')
      % Inclined paths
      Re = 8500; % [km] effective Earth radius, ITU_R P.834
      phi1 = phi; % elevation angle at h1
      phi2 = arccos((Re + h2)/(Re + h2)*cos(phi1));
      x1o = tan(phi1)*sqrt((Re+h1)/h_o);   x2o = tan(phi2)*sqrt((Re+h2)/ho);
      x1w = tan(phi1)*sqrt((Re+h1)/h_w);   x2w = tan(phi2)*sqrt((Re+h1)/hw);
      F = @(x) 1 / (0.661*x + 0.339*sqrt(x^2+5.51));
      
      A = gam_o*sqrt(h_o)*(sqrt(Re+h1)*F(x1o)*exp(-h1/h_o)/cos(phi1) - sqrt(Re+h2)*F(x2o)*exp(-h2/ho)/cos(phi2)  ) + ...
          gam_w*sqrt(h_w)*(sqrt(Re+h1)*F(x1w)*exp(-h1/h_w)/cos(phi1) - sqrt(Re+h2)*F(x2w)*exp(-h2/hw)/cos(phi2)  ); %[dB]
       end
end


end

