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

clc;
close all;
clear all;

set(0,'defaultTextInterpreter','latex')
fontSize = 20;
set(0,'DefaultAxesFontSize', fontSize);

addpath('../savedData')
addpath('../functions')

load ../savedData/attenuationMap.mat

run('System')
run('GroundStation.m')

a = 7.204090326167004e3;     % Semi-major axis [km] 
R_E = 6378.137;              % Earth's equatorial radius [km]

% Maximum trasmitter-receiver distance [m] (i.e. most critical condition at 5°)
d_max = 1e+3*( - R_E*sin(th_el) + sqrt((R_E*sin(th_el))^2 + (a-R_E)^2 + 2*R_E*(a-R_E)) );

% Transmitter
TX.P = 10 ^ ((27-30)/10);   % [W] Trasmitted power (27 dBm RF power output)
TX.G = db2pow(6);           % Gain of TX antenna [W]
TX.L = db2pow(-3);          % Additional losses in the transmitter chain [W]

% Receiver
RX.P_acc = 5e-3;            % Pointing  [deg]
RX.L = db2pow(-3);          % Additional losses in the receiver chain [W]
Eb_N0_min = 7;     


c = 3*1e+08;                % Speed of light [m/s]
lambda = c/(f*1e+9);        % Wavelength [m]
k_B = 1.379*1e-23;          % Boltzmann constant [J/K]


RX.Lp = db2pow(- 12*(GS1.D *RX.P_acc / (70*lambda))^2) ; %  RX pointing losses [dB]
TX.Lp = db2pow(-0.1); %  TX pointing losses [dB]

A_atm = db2pow(- A_tot_GS(1,1)) ; % Atmospheric attenuation at 5 deg, critical angle

L_FP =  (4*pi*d_max*f*1e9/c)^-2; % Free path loss

% Eb/N0 
Eb_N0 = TX.P * TX.G * TX.Lp * TX.L * db2pow(GS1.G_T) * RX.Lp * RX.L * A_atm * L_FP / (k_B*(R_D*1e+9));

Eb_N0 = pow2db(Eb_N0);

link_margin = Eb_N0 - Eb_N0_min;





