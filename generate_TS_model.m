function [Ac,Bc] = generate_TS_model(params,Ts,h_vals)
%GENERATE_TS_MODEL Generate a discrete-time Takagi–Sugeno (TS) model of a two-tank system.
%
%   [Ac,Bc] = generate_TS_model(params,Ts,h_vals)
%
%   The local continuous-time models are obtained by Jacobian linearization
%   of the nonlinear Torricelli outflow dynamics around the operating points
%   (h1_i, h2_j) and then discretized using ZOH.
%
%   Assumptions:
%     - Inputs are PWM duty cycles u1,u2 in [0,1]
%     - Inflows: Qin1 = Q_max,1*u1, Qin2 = Q_max,2*u2
%       (this implementation uses params.Q_max for both channels; if pumps differ,
%       provide params.Q_max as a 2-element vector and adjust B_ij accordingly)
%
%   Inputs:
%     params : struct with fields including physical parameters such as At1, At2, k1, k2, Q_max
%     Ts     : sampling time [s]
%     h_vals : grid of operating points for h1 and h2 (positive, avoid 0)
%
%   Outputs:
%     Ac, Bc : cell arrays (N^2 x 1) of discrete-time local models
%
%   Note:
%     Discretization is performed without ss/c2d to avoid Control System Toolbox
%     dependency (ZOH via augmented matrix exponential).

%% Basic parameters
H1 = params.H1;
H2 = params.H2;
D1 = params.D1;
D2 = params.D2;
g = params.g;
R1 = params.R1;
R2 = params.R2;
At1 = params.At1;
At2 = params.At2;
Q_max = params.Q_max;
d_out1 = params.d_out1;
d_out2 = params.d_out2;
A_out1 = params.A_out1;
C_out1 = params.C_out1;
A_out2 = params.A_out2;
C_out2 = params.C_out2;
k1 = params.k1;
k2 = params.k2;

N = length(h_vals);

% Container for models
Ac=cell(N*N,1);
Bc=cell(N*N,1);

% iterating variable
model_idx = 1;

for i = 1:N
    for j = 1:N
        h1 = h_vals(i);
        h2 = h_vals(j);
        sqrt_h1 = sqrt(h1);
        sqrt_h2 = sqrt(h2);

        % Based on the Jacobian of the nonlinear model
        % Contineous-time model linearized around (h1, h2)
        A_ij = [-0.5*k1 / (At1 * sqrt_h1),       0;
            0.5 * k1 / (1*At2 * sqrt_h1),  -0.5 * k2 / (1*At2 * sqrt_h2)];

        B_ij = [Q_max / At1,      0;
            0,         Q_max / At2];

        % model discretization
        % ZOH discretization
        n = size(A_ij,1);
        r = size(B_ij,2);
        % Van Loan / augmented matrix exponential
        M  = [A_ij, B_ij;
              zeros(r,n), zeros(r,r)];
        Md = expm(M*Ts);
        Ad = Md(1:n, 1:n);
        Bd = Md(1:n, n+1:n+r);

        % save as cell
        Ac{model_idx,1} = Ad;
        Bc{model_idx,1} = Bd;

        model_idx = model_idx + 1;
    end
end

%% Additionally save it for future useage
% save('TSmodel.mat', 'Ac', 'Bc', 'h_vals', 'Ts');
%
% disp('✅ Takagi-Sugeno model saved as: TSmodel.mat');

%%
end