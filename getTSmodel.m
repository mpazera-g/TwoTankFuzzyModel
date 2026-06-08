function [Ad,Bd,mu_vec] = getTSmodel(h1,h2,Ac,Bc,h_vlas,type,sigma)
%GETTSMODEL interpolates Takagi–Sugeno model matrices using specified membership functions
%
%   [Ad,Bd,mu_vec] = getTSmodel(h1,h2,A_cell,B_cell,grid_vals,type,sigma)
%
% Inputs:
%   h1, h2      - current values of premise variables
%   Ac          - cell array with local A matrices (size N^2 x 1)
%   Bc          - cell array with local B matrices (size N^2 x 1)
%   h_vlas      - vector of premise points (e.g. [0.2, 0.5, 1.0])
%   type        - (optional) 'triangular' [default] or 'gaussian'
%   sigma       - (required if type = 'gaussian')
%
% Outputs:
%   Ad, Bd      - interpolated system matrices
%   mu_vec      - activation weights (N^2 x 1)

    if nargin < 6
        type = 'triangular';  % default type
    end

    N = length(h_vlas);          % number of rules per variable
    num_rules = N * N;

    %% --- Compute membership weights ---
    switch lower(type)
        case 'triangular'
            mu_h1 = interpMembership(h1, h_vlas);  % 1 x N
            mu_h2 = interpMembership(h2, h_vlas);  % 1 x N
        case 'gaussian'
            if nargin < 7
                error('For gaussian membership, provide sigma.');
            end
            mu_h1 = exp(-((h1 - h_vlas).^2) / (2 * sigma^2));
            mu_h2 = exp(-((h2 - h_vlas).^2) / (2 * sigma^2));
        otherwise
            error('Unknown membership function type: %s', type);
    end

    % Outer product to form rule activations
    mu_mat = mu_h2' * mu_h1;          % size N x N
    mu_vec = reshape(mu_mat, [], 1);  % size N^2 x 1

    % Normalize weights
    mu_vec = mu_vec / (sum(mu_vec) + eps);  % avoid division by zero

    %% --- Interpolate system matrices ---
    Ad = zeros(size(Ac{1}));
    Bd = zeros(size(Bc{1}));

    for i = 1:num_rules
        Ad = Ad + mu_vec(i) * Ac{i};
        Bd = Bd + mu_vec(i) * Bc{i};
    end
end

%% --- Local helper: strictly interpolative triangular membership ---
function mu = interpMembership(val, grid)
% Returns vector of memberships with max 2 non-zero values (linear interpolation)
    n = length(grid);
    mu = zeros(1, n);

    % Below range
    if val <= grid(1)
        mu(1) = 1;
        return;
    end

    % Above range
    if val >= grid(end)
        mu(end) = 1;
        return;
    end
    
    % Inside range
    idx = find(grid <= val, 1, 'last');
    if grid(idx) == val
        mu(idx) = 1;
    else
        a = grid(idx);
        b = grid(idx + 1);
        mu(idx)     = (b - val) / (b - a);
        mu(idx + 1) = (val - a) / (b - a);
    end
end
