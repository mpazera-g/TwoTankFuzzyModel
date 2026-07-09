function [Ad,Bd,cd] = generate_TS_model(params,Ts,h_vals)
%GENERATE_TS_MODEL Generate a discrete-time affine Takagi-Sugeno (TS) model
%of a nonlinear two-tank system.
%
%   [Ad,Bd] = generate_TS_model(params,Ts,h_vals)
%   [Ad,Bd,cd] = generate_TS_model(params,Ts,h_vals)
%
%   The local continuous-time models are obtained by Jacobian linearization
%   of the nonlinear Torricelli outflow dynamics around the operating points
%   (h1_i,h2_j). When the local model is written in absolute variables, it
%   has the affine form
%
%       xdot = A_ij*x + B_ij*u + c_ij.
%
%   Each local affine model is discretized using exact ZOH via an augmented
%   matrix exponential, resulting in
%
%       x(k+1) = Ad_ij*x(k) + Bd_ij*u(k) + cd_ij.
%
%   Backward compatibility:
%     - If the function is called with two outputs, only Ad and Bd are
%       returned, as in the original implementation.
%     - If the function is called with three outputs, cd is also returned.
%
%   Inputs:
%     params : struct with fields At1, At2, k1, k2, Q_max, etc.
%     Ts     : sampling time [s]
%     h_vals : grid of positive operating points for h1 and h2
%
%   Outputs:
%     Ad, Bd : cell arrays (N^2 x 1) of discrete-time local matrices
%     cd     : cell array  (N^2 x 1) of discrete-time affine terms
%
%   Note:
%     h_vals must contain strictly positive values, because the Jacobian of
%     sqrt(h) is singular at h = 0.

    %% Input checks
    if nargin < 3
        error('generate_TS_model requires params, Ts and h_vals.');
    end

    if Ts <= 0
        error('Sampling time Ts must be positive.');
    end

    if any(h_vals <= 0)
        error('All operating points in h_vals must be strictly positive.');
    end

    %% Physical parameters
    At1 = params.At1;
    At2 = params.At2;
    k1  = params.k1;
    k2  = params.k2;

    % Pump gains. A scalar Q_max means identical pumps. A two-element vector
    % allows different maximum flow rates for the two pumps.
    Q_max = params.Q_max;
    if isscalar(Q_max)
        Q_max1 = Q_max;
        Q_max2 = Q_max;
    elseif numel(Q_max) == 2
        Q_max1 = Q_max(1);
        Q_max2 = Q_max(2);
    else
        error('params.Q_max must be either a scalar or a two-element vector.');
    end

    %% Containers for local models
    N = numel(h_vals);
    num_rules = N*N;

    Ad = cell(num_rules,1);
    Bd = cell(num_rules,1);
    cd = cell(num_rules,1);

    %% Generate local affine models and discretize them
    model_idx = 1;

    for i = 1:N
        for j = 1:N
            h1 = h_vals(i);
            h2 = h_vals(j);

            sqrt_h1 = sqrt(h1);
            sqrt_h2 = sqrt(h2);

            % Continuous-time Jacobian matrix A_ij
            A_ij = [ -0.5*k1/(At1*sqrt_h1),                  0;
                      0.5*k1/(At2*sqrt_h1), -0.5*k2/(At2*sqrt_h2) ];

            % Continuous-time input matrix B_ij. It is constant here because
            % pump inflows enter the model linearly as Qin,j = Q_max,j*u_j.
            B_ij = [ Q_max1/At1,          0;
                              0, Q_max2/At2 ];

            % Continuous-time affine term c_ij for absolute variables.
            % The nonlinear vector field is f(x,u) = g(x) + B*u, hence
            % c_ij = f(x_ij,u_ij) - A_ij*x_ij - B_ij*u_ij = g(x_ij)-A_ij*x_ij.
            g_ij = [ -k1*sqrt_h1/At1;
                     (k1*sqrt_h1 - k2*sqrt_h2)/At2 ];

            x_ij = [h1; h2];
            c_ij = g_ij - A_ij*x_ij;

            % Exact ZOH discretization of the affine continuous-time system:
            %   xdot = A*x + B*u + c
            % by augmenting B with the constant input channel c.
            n = size(A_ij,1);
            r = size(B_ij,2);

            M = [A_ij, B_ij, c_ij;
                 zeros(r+1,n+r+1)];

            Md = expm(M*Ts);

            Ad{model_idx,1} = Md(1:n,1:n);
            Bd{model_idx,1} = Md(1:n,n+1:n+r);
            cd{model_idx,1} = Md(1:n,n+r+1);

            model_idx = model_idx + 1;
        end
    end
end
