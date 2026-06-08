%main.m
% Nonlinear simulation of the two-tank benchmark
% based on the dynamical equations with:
%   - TS model generation (local linear models + ZOH discretization)
%   - TS interpolation (triangular or gaussian membership)
%   - LMI-based (common) state-feedback controller design
%   - LMI-based (common) state observer design
%
% Notes:
%   - Inputs are PWM duty cycles u in [0,1].
%   - Optional inter-tank transport delay is implemented in the simulator
%     via a finite buffer on Q_out,1 -> Q_in,2 (not included in TS matrices).

%%
clear; close all; clc;

%% Simulation parameters
Ts=0.1;                 % [s] sampling time
T_end=600;              % [s] simulation time
nt=round(T_end/Ts) + 1; % no. of steps
t=(0:nt-1)*Ts;          % time vector

%% get two-Tank system data, its basic physical parameters of the two-tank system
T=getTanks();

%% grid of fuzzy values h1, h2
h_low=0.1;
h_high=0.9;
noFuzzySets=5;
h_vals=linspace(h_low,h_high,noFuzzySets);  % Caution: Please do not use 0 as the lower value to avoid dividing by 0

%% TS model generation
[Ad,Bd]=generate_TS_model(T,Ts,h_vals);

%% state space output matrix
C=eye(2);   % output matrix of the system like y(k)=C*x(k); 2 state variable are measurable

%% size of the matrices
n=size(Ad{1},1);        % A=nxn
r=size(Bd{1},2);        % B=nxr
m=size(C,1);            % C=mxn
nts=size(Ad,1);         % no. of TS model matrices

%% check controlability and observability
[co,ob]=checkContObs(Ad,Bd,C);

%% state-feedback LMI-based controller design
K=designController(Ad,Bd);
Cc=[1,0; 0,1]; % controlled output

%% LMI-based state observer design
L=designObserver(Ad,C);

%% reference signals to which the system should by guided
ref=zeros(r,nt);
ref(2,:)=0.5;                   % constant reference for the lower tank
omega=2*pi/nt/Ts;               % one sine period over the simulation time
ref(1,:)=0.3*sin(omega*t)+0.5;  % sine-like reference for the upper tank

%% container variables for output TS matrices
A_sim=cell(1,nt-1);
B_sim=cell(1,nt-1);
w=zeros(nts, nt);
A_est=cell(1,nt-1);
B_est=cell(1,nt-1);
w_est=zeros(nts,nt);

%%
T_d=3.0;                                % water Q_out,1-> Q_in,2 delay [s]
delay_steps=round(T_d/Ts);              % no. of delay steps
Q_out1_hist=zeros(1, delay_steps+1);    % delay buffer

%% container for simulation vectors
x=zeros(n,nt);            % state of the system
x(:,1)=[0;0];             % initial condition
h1=zeros(1, nt);          % water level in the first tank
h2=zeros(1, nt);          % water level in the first tank
h1(1)=0;                  % initial conditions
h2(1)=0;
y=zeros(m,nt);            % output of the system
u1=ones(1, nt)*1;         % PWM for the first pomp (0–1)
u2=ones(1, nt)*1;         % PWM for the second pomp (0–1)
u=[u1;u2];
xest=zeros(n,nt);         % state estimate
xest(:,1)=[0.1;0.05];     % state estimate initial condition
Q_out1_log=zeros(1, nt);  % whole outflow from the upper tank log
Q_out2_log=zeros(1, nt);  % whole outflow from the lower tank log
Q_in2_log=zeros(1, nt);   % delayed inflow to the lower tank log

%% simulation-main loop
k=1;
while k < nt
    % get current water level
    h1_k=x(1,k);
    h2_k=x(2,k);

    % collision avoidance
    sqrt_h1=sqrt(max(h1_k, 0));
    sqrt_h2=sqrt(max(h2_k, 0));

    sigma=0.25;
    % [A_sim{k},B_sim{k},w(:,k)]=getTSmodel(h1_k,h2_k,Ac,Bc,h_vals,'gaussian',sigma); % get current A and B matrices based on TS model
    [A_sim{k},B_sim{k},w(:,k)]=getTSmodel(h1_k,h2_k,Ad,Bd,h_vals); % get current A and B matrices based on TS model

    M  = -Cc/(A_sim{k}-B_sim{k}*K-eye(n))*B_sim{k}; % auxilary variable for pre-filter Kr
    Kr = M \ eye(size(M,1)); % calculate a pre-filter matrix to the reference vector
    u(:,k)=-K*x(:,k)+Kr*ref(:,k); % calculate control input
    % u(:,k)=[1;1];   % open-loop control
    u(:,k) = min(max(u(:,k),0),1); % control saturation in [0,1]

    % calculating the water levels with using nonlinear model
    dh1=(T.Q_max*u(1,k)-T.k1*sqrt_h1)/T.At1;
    Q_out1_now=T.k1*sqrt_h1; % calculate Q_out1
    Q_out1_hist=[Q_out1_now, Q_out1_hist(1:end-1)]; % update the buffer
    Q_in2_delayed=Q_out1_hist(end);    % get delayed flow
    dh2=(Q_in2_delayed+T.Q_max*u(2,k)-T.k2*sqrt_h2)/T.At2; % now update dh2 with delayed flow

    % Euler discretization
    h1(k+1)=h1_k+Ts*dh1;
    h2(k+1)=h2_k+Ts*dh2;
    x(:,k+1)=[h1(k+1);h2(k+1)];
    y(:,k)=C*x(:,k);

    [A_est{k},B_est{k},w_est(:,k)]=getTSmodel(xest(1,k),xest(2,k),Ad,Bd,h_vals); % get current A and B matrices based on TS model
    xest(:,k+1)=A_est{k}*xest(:,k)+B_est{k}*u(:,k)+L*(y(:,k)-C*xest(:,k)); % state estimation

    k=k+1; % update time instance
end

%% Plot results: water levels and estimates
figure;
for i=1:n
    subplot(n,1,i);
    hold on;
    plot(t, x(i,:), 'b', 'LineWidth', 2);
    plot(t, ref(i,:), '--r', 'LineWidth', 2);
    plot(t, xest(i,:), '--g', 'LineWidth', 2);
    grid on;
    axis([0 T_end 0 max(x(i,:))+0.1]);
    ylabel(strcat('h_',num2str(i),' [m]'));
    xlabel('Time [s]');
    legend('water level', 'reference','state estimate','Location','Best');
end

%% Plot results: control input
figure;
for i=1:r
    subplot(r,1,i);
    plot(t, u(i,:), 'b', 'LineWidth', 1.5); hold on;
    xlabel('Time [s]');
    ylabel(strcat('u_',num2str(i)));
    axis([0 T_end 0 1.1]);
    grid on;
end

%% Plot results: reference
figure;
for i=1:n
    subplot(n,1,i);
    plot(t, ref(i,:), 'b', 'LineWidth', 1.5); hold on;
    xlabel('Time [s]');
    ylabel(strcat('r_',num2str(i)));
    axis([0 T_end 0 1.1]);
    grid on;
end

%%
mus=zeros(1,nt-1);
for i=1:nt-1
    mus(1,i)=sum(w(:,i));
end

%% Plot results: plotmembership functions
h_plot = linspace(0, 1, 1000);   % y axis
figure;
% --- triangular MF ---
subplot(2,1,1);
hold on;
title('Triangular Membership Functions');
xlabel('water level','Interpreter','Latex'); ylabel('$\mu$','Interpreter','Latex');
for i = 1:length(h_vals)
    mu = zeros(size(h_plot));
    if i == 1
        % Left edge, then falling
        a = h_vals(1);
        b = h_vals(2);
        mu(h_plot <= a) = 1;
        idx = h_plot > a & h_plot <= b;
        mu(idx) = (b - h_plot(idx)) / (b - a);
    elseif i == length(h_vals)
        % Rising, then right edge
        a = h_vals(end - 1);
        b = h_vals(end);
        idx = h_plot >= a & h_plot < b;
        mu(idx) = (h_plot(idx) - a) / (b - a);
        mu(h_plot >= b) = 1;
    else
        % Classical triangular MF
        a = h_vals(i - 1);
        b = h_vals(i);
        c = h_vals(i + 1);
        idx1 = h_plot >= a & h_plot <= b;
        idx2 = h_plot >= b & h_plot <= c;
        mu(idx1) = (h_plot(idx1) - a) / (b - a);
        mu(idx2) = (c - h_plot(idx2)) / (c - b);
    end
    plot(h_plot, mu, 'LineWidth', 2);
end
ylim([0 1.1]); grid on;


% --- Gaussian MF ---
subplot(2,1,2);
hold on;
title('Gaussian Membership Functions');
xlabel('water level','Interpreter','Latex'); ylabel('$\mu$','Interpreter','Latex');
for i = 1:length(h_vals)
    mu = exp(-((h_plot - h_vals(i)).^2) / (2 * sigma^2));
    plot(h_plot, mu, 'LineWidth', 2);
end
ylim([0 1.1]); grid on;
legend(arrayfun(@(x) sprintf('h = %.2f', x), h_vals, 'UniformOutput', false));

%%

%%
