function [gam,h] = gas_specAtt(f,gas,GS,Map,Table)

rho = interp2(Map.LON_rho, Map.LAT_rho, Map.rho1_map, GS.lon, GS.lat);
[T,p]=read_mprof(GS.lat,GS.lon);
theta = 300/T;
e = rho*T / 216.7;
rp = (p+e)/(1013.25);

% Nd is the dry continuum due to pressure-induced nitrogen absorption and
% the Debye spectrum (ITU-R P.676-8 pag.8)
d = 5.6*1e-04*(p+e)*theta^0.8; % width parameter for the Debye spectrum
Nd = f*p*theta^2*((6.14*1e-05/(d*(1+(f/d)^2)))+...
    (1.4*1e-12*p*theta^1.5/(1+1.9*1e-05*f^1.5)));

% OXYGEN
if strcmp(gas,'oxygen')
    fo=Table.ox(1:end,1);
    a1=Table.ox(1:end,2);
    a2=Table.ox(1:end,3);
    a3=Table.ox(1:end,4);
    a4=Table.ox(1:end,5);
    a5=Table.ox(1:end,6);
    a6=Table.ox(1:end,7);
    N = length(fo);
    
    N_o_vec = [];
    % ITU-R P.676-8
    for i = 1:N
        % Correction factor whihc arises due to interference effects in
        % Oxygen lines
        delta_o = (a5(i) + a6(i)*theta)*1e-04*(p+e)*theta^0.8;
        % Line width
        Df_o = a3(i)*1e-04*(p*theta^(0.8-a4(i)) + 1.1*e*theta);
        Df_o = sqrt(Df_o^2 + 2.25*1e-06);
        % Line shape factors
        F_o = f/fo(i)*( ((Df_o-delta_o*(fo(i)-f))/((fo(i)-f)^2 + Df_o^2)) + ...
            ((Df_o-delta_o*(fo(i)+f))/((fo(i)+f)^2 + Df_o^2)) );
        % Strength of the i-th oxygen or water vapour line
        S_o = a1(i)*1e-07*p*theta^3 * exp(a2(i)*(1-theta));
        % Imaginary parts of the frequency-dependent complex refractivities
        N_o_vec = [N_o_vec;S_o*F_o];
    end
    % Extend over all the lines
    N_o = sum(N_o_vec)+ Nd;
    gam = 0.1820*f*N_o; % [dB/km] Specific attenuations due to dry air
    
    % h_o: equivalent height attributable to the oxygen component of gaseous
    % attenuation. constraint: h_o<10.7*rp^0.3, f<70 GHz
    h_o_lim = 10.7*rp^0.3;
    A = 0.7832 + 0.00709*(T-273.15); %(34)
    t1 = 5.1040/(1 + 0.066*rp^(-2.3)) * exp(-((f-59.7)/(2.87 + 12.4*exp(-7.9*rp)))^2); %(31)
    
    % values for f and c Table 3 pag.25 
    c_vec = [0.1597 0.1066 0.1325 0.1242 0.0938 0.1448 0.1374];
    f_vec = [118.750334 368.498246 424.763020 487.249273 715.392902 773.839490 834.145546]; % [GHz]
    
    t2= sum( c_vec.*exp(2.12.*rp)./((f-f_vec).^2 + 0.025*exp(2.2*rp))) ;
    
    t3 = 0.0114*f*(15.02*f^2 - 1353*f + 5.333e+04)/...
        ((1+0.14*rp^(-2.6))*(f^3-151.3*f^2+9629*f-6803));
    % h0 : equivalent height attributable to the oxygen component of
    % gaseous attenutation
    h0 = 6.1*A*(1+t1+t2+t3) / (1 + 0.17*rp^(-1.1));
    h = h0;

% WATER VAPOUR    
elseif strcmp(gas,'waterVapour')
    
    fw=Table.wv(1:end,1);
    b1=Table.wv(1:end,2);
    b2=Table.wv(1:end,3);
    b3=Table.wv(1:end,4);
    b4=Table.wv(1:end,5);
    b5=Table.wv(1:end,6);
    b6=Table.wv(1:end,7);
    M = length(fw);
    
   
    N_w_vec = [];
    for i = 1:M
        delta_w = 0; %(7)
        Df_w = b3(i)*1e-04*(p*theta^(b4(i)) + b5(i)*e*theta^b6(i)); %(6a)
        Df_w = 0.535*Df_w+ sqrt(0.217*Df_w^2 + (2.1316*1e-12*fw(i)^2)/theta); %(6b)
        F_w = f/fw(i)*( ((Df_w-delta_w*(fw(i)-f))/((fw(i)-f)^2 + Df_w^2)) + ...
            ((Df_w-delta_w*(fw(i)+f))/((fw(i)+f)^2 + Df_w^2)) ); %(5)
        S_w = b1(i)*1e-01*e*theta^3.5 * exp(b2(i)*(1-theta)); %(3)
        N_w_vec = [N_w_vec;S_w*F_w];
    end
    N_w = sum(N_w_vec);
    
  
    gam = 0.1820*f*N_w ;% [dB/km]specific attenuations due to water vapour
        
    % h_w: equivalent height attributable to the water vapour component of
    % gaseous attenuation
    A = 1.9298 - 0.04166*(T-273.15)+ 0.0517*rho;
    B = 1.1674 - 0.00622*(T-273.15)+ 0.0063*rho;
    sigmaw = (1.013)/(1+exp(-8.6*(rp-0.57)));
    
    fi=Table.wv_height(1:end,1);
    a=Table.wv_height(1:end,2);
    b=Table.wv_height(1:end,3);
    
    hw = A + B*sum( a.*sigmaw./((f-fi).^2 + b.*sigmaw));
    h = hw;
     
else
    error('Invalid INPUT gas, type "oxygen" or "waterVapour" ')
end


end

