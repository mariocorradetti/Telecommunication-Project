function [T,p]=read_mprof(lat,lon)
% extracts the mean monthly vertical profiles at 00, 06, 12 and 18 UTC
% of air temperature, total pressure 
% of the pixel located at latitude and longitude indexes ilat and ilon
% 
% It uses input data derived by ESA from analysis of ECMWF ERA15 product.
%
% INPUT
%	ilat= latitude index,    range (1..121) 1=90 2=88.5 121=-90 [deg] 
%   ilong = longitude index, range=(1..240) 1=0 2=1.5 240=358.5 [deg]
% 	N.B. 
%       ilat = ilong = -1 : close the input files
%
% OUTPUT
%  tk: Temperature of air [K],                   size=nlev X nhours*nmonths
%  pr: Air total Pressure [hPa],                 size=nlev X nhours*nmonths
%
% N.B.
% nlev=32 number of vertical levels
% nhours=4 number of hours 1=00, 2=06, 3=12 and 4=18 UTC 
% nmonths= 12 number of months 1=Jan 12=Dec
%
% e.g. tk(:,1:4),hm(:,1) is the mean vertical profile of 
%         air temperature at 00, 06, 12 and 18 UTC in January 
ilat = floor((90-lat)/1.5)+1;
ilong = floor(lon/1.5)+1;
% Relase 1.0 11/2003
% by Antonio Martellucci, Giulio Blarzino, ESA/ESTEC
    tk=read_tk(ilat,ilong);
	pr=read_pr(ilat,ilong);
    
    T = mean(tk(1,1:end));
    p = mean(pr(1,1:end));
    


return
