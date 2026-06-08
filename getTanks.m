function T=getTanks()
%GETTANKS Parameters of the benchmark two-tank system.
%
%   T = getTanks() returns a struct with physical and hydraulic parameters used
%   in the nonlinear simulator and TS model generation.

%% create container for data
T=struct;
%% put data into struct
%% physical consts
T.H1=2;             % upper tank height [m]
T.D1=0.2;           % upper tank diameter [m]
T.H2=2;             % lower tank height [m]
T.D2=0.2;           % lower tank diameter [m]
T.g=9.81;           % gravity [m/s^2]
%% Cross-sectional are of the tank
T.R1=T.D1 / 2;            % radius of the upper tank [m]
T.R2=T.D2 / 2;            % radius of the lower tank [m]
T.At1=pi*T.R1^2;        % cross-sectional area of the upper tank [m^2]
T.At2=pi*T.R2^2;        % cross-sectional area of the lower tank [m^2]
%% pumps parameters
T.Q_max=1e-3;           % pump efficiency [m^3/s]; 1e-3[m^3/s] ~ 60 [l/min]
%% outflow parameters
T.d_out1=0.02;                      % opening hole diameter [m] of the upper tank
T.A_out1=pi*(T.d_out1/2)^2;         % cross-sectional area of the opening hole [m^2] of the upper tank
T.C_out1=0.72;                      % outflow coefficient of the upper tank
T.d_out2=0.023;                     % opening hole diameter [m] of the lower tank
T.A_out2=pi*(T.d_out2/2)^2;         % cross-sectional area of the opening hole [m^2] of the lower tank
T.C_out2=1;                         % outflow coefficient of the lower tank
T.k1=T.C_out1*T.A_out1*sqrt(2*T.g); % outflow constant [m^3/s/sqrt(m)] of the upper tank
T.k2=T.C_out2*T.A_out2*sqrt(2*T.g); % outflow constant [m^3/s/sqrt(m)] of the lower tank
end