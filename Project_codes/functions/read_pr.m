function [pr]=read_pr(ilat,ilong)
% extracts the mean monthly vertical profiles at 00, 06, 12 and 18 UTC
% of air  total pressure a 
% of the pixel located at latitude and longitude indexes ilat and ilon
% 
% It uses input data derived by ESA from analysis of ECMWF ERA15 product.
%
% INPUT
%	ilat= latitude index,    range (1..121) 1=90 2=88.5 121=-90 [deg] 
%   ilong = longitude index, range=(1..240) 1=0 2=1.5 240=358.5 [deg]
% 	N.B. 
%       ilong=241=1  
%       ilat = ilong = -1 : close the input files
%
% OUTPUT
%  pr: Air total Pressure [hPa],                 size=nlev X nhours*nmonths
%
% N.B.
% nlev=32 number of vertical levels
% nhours=4 number of hours 1=00, 2=06, 3=12 and 4=18 UTC 
% nmonths= 12 number of months 1=Jan 12=Dec
%
% e.g. pr(:,1:4),hm(:,1) is the mean air total pressure profile 
%      at 00, 06, 12 and 18 UTC in January 

% Relase 1.0 11/2003
% by Antonio Martellucci, Giulio Blarzino, ESA/ESTEC

persistent fid 

if ((ilat == -1) & (ilong == -1))
   if ~isempty(fid)
      fclose(fid(1));
      fclose(fid(2));
      fclose(fid(3));
      fclose(fid(4));
      clear fid;
   end
   pr=[];
else
   if isempty(fid)
		fid(1)=fopen('atm_properties/pres_00.bin','r','ieee-be');
		fid(2)=fopen('atm_properties/pres_06.bin','r','ieee-be');
		fid(3)=fopen('atm_properties/pres_12.bin','r','ieee-be');
		fid(4)=fopen('atm_properties/pres_18.bin','r','ieee-be');
	end   
	x=zeros(48,32);
	lrec=1536;
	ipos=((ilong-1)+(ilat-1)*241).*lrec;

	fseek(fid(1),ipos,'bof');
	x(1:4:48,:)=fread(fid(1),[12,32],'single');
	fseek(fid(2),ipos,'bof');
	x(2:4:48,:)=fread(fid(2),[12,32],'single');
	fseek(fid(3),ipos,'bof');
	x(3:4:48,:)=fread(fid(3),[12,32],'single');
	fseek(fid(4),ipos,'bof');
	x(4:4:48,:)=fread(fid(4),[12,32],'single');
	pr=x';
end

return
