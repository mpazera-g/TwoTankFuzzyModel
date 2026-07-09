function [Ad,Bd,varargout] = getTSmodel(h1,h2,Ac,Bc,varargin)
%GETTSMODEL Interpolate Takagi-Sugeno model matrices using membership functions.
%
%   Backward-compatible syntax without affine terms:
%     [A_k,B_k,mu_vec] = getTSmodel(h1,h2,Ac,Bc,h_vals)
%     [A_k,B_k,mu_vec] = getTSmodel(h1,h2,Ac,Bc,h_vals,type,sigma)
%
%   New syntax with affine terms:
%     [A_k,B_k,c_k,mu_vec] = getTSmodel(h1,h2,Ac,Bc,cc,h_vals)
%     [A_k,B_k,c_k,mu_vec] = getTSmodel(h1,h2,Ac,Bc,cc,h_vals,type,sigma)
%
%   Inputs:
%     h1, h2 : current values of premise variables
%     Ac     : cell array with local A matrices (N^2 x 1)
%     Bc     : cell array with local B matrices (N^2 x 1)
%     cc     : optional cell array with local affine terms (N^2 x 1)
%     h_vals : vector of positive premise grid points
%     type   : optional membership type: 'triangular' [default] or 'gaussian'
%     sigma  : Gaussian width parameter; required for 'gaussian'
%
%   Outputs:
%     A_k, B_k : interpolated discrete-time matrices
%     c_k      : interpolated affine term, returned only when four outputs
%                are requested
%     mu_vec   : normalized activation weights (N^2 x 1)
%
%   Compatibility rule:
%     - With three outputs, the third output is mu_vec, exactly as in the
%       original implementation.
%     - With four outputs, the third output is c_k and the fourth output is
%       mu_vec.

    %% Parse inputs
    if nargin < 5
        error('Not enough input arguments.');
    end

    firstOptional = varargin{1};

    if iscell(firstOptional)
        % New syntax: getTSmodel(h1,h2,Ac,Bc,cc,h_vals,...)
        cc = firstOptional;

        if numel(varargin) < 2
            error('When affine terms are provided, h_vals must be the next argument.');
        end

        h_vals = varargin{2};
        remaining = varargin(3:end);
    else
        % Old syntax: getTSmodel(h1,h2,Ac,Bc,h_vals,...)
        h_vals = firstOptional;
        remaining = varargin(2:end);

        % No affine terms supplied: use zero affine terms.
        cc = cell(size(Ac));
        for i = 1:numel(Ac)
            cc{i} = zeros(size(Ac{i},1),1);
        end
    end

    if isempty(remaining)
        type = 'triangular';
    else
        type = remaining{1};
    end

    if numel(remaining) >= 2
        sigma = remaining{2};
    else
        sigma = [];
    end

    %% Basic checks
    N = numel(h_vals);
    num_rules = N*N;

    if numel(Ac) ~= num_rules || numel(Bc) ~= num_rules || numel(cc) ~= num_rules
        error('The numbers of A, B and c local models must all be equal to numel(h_vals)^2.');
    end

    %% Compute membership weights
    switch lower(type)
        case 'triangular'
            mu_h1 = interpMembershipTriangular(h1,h_vals);
            mu_h2 = interpMembershipTriangular(h2,h_vals);

        case 'gaussian'
            if isempty(sigma) || sigma <= 0
                error('For gaussian membership, provide a positive sigma.');
            end

            mu_h1 = exp(-((h1 - h_vals).^2)/(2*sigma^2));
            mu_h2 = exp(-((h2 - h_vals).^2)/(2*sigma^2));

        otherwise
            error('Unknown membership function type: %s', type);
    end

    % Rule ordering is consistent with generate_TS_model.m:
    % for i = 1:N      % h1 index
    %   for j = 1:N    % h2 index
    %      idx = idx+1
    %   end
    % end
    mu_vec = zeros(num_rules,1);
    idx = 1;
    for i = 1:N
        for j = 1:N
            mu_vec(idx) = mu_h1(i)*mu_h2(j);
            idx = idx + 1;
        end
    end

    weightSum = sum(mu_vec);
    if weightSum <= eps
        warning('All membership weights are numerically zero. Uniform weights are used.');
        mu_vec = ones(num_rules,1)/num_rules;
    else
        mu_vec = mu_vec/weightSum;
    end

    %% Interpolate A, B and c
    Ad = zeros(size(Ac{1}));
    Bd = zeros(size(Bc{1}));
    cd = zeros(size(cc{1}));

    for i = 1:num_rules
        Ad = Ad + mu_vec(i)*Ac{i};
        Bd = Bd + mu_vec(i)*Bc{i};
        cd = cd + mu_vec(i)*cc{i};
    end

    %% Return outputs with backward compatibility
    if nargout <= 3
        % Original behavior: [A_k,B_k,mu_vec]
        varargout{1} = mu_vec;
    else
        % New behavior: [A_k,B_k,c_k,mu_vec]
        varargout{1} = cd;
        varargout{2} = mu_vec;
    end
end

%% Local helper: boundary-saturated triangular membership
function mu = interpMembershipTriangular(val,grid)
%INTERPMEMBERSHIPTRIANGULAR Return triangular membership degrees.
%
% The function returns at most two non-zero values inside the grid. Outside
% the grid, boundary saturation is used:
%   val <= grid(1)   -> first membership equals 1
%   val >= grid(end) -> last membership equals 1

    n = numel(grid);
    mu = zeros(1,n);

    if val <= grid(1)
        mu(1) = 1;
        return;
    end

    if val >= grid(end)
        mu(end) = 1;
        return;
    end

    idx = find(grid <= val,1,'last');

    if abs(grid(idx)-val) < 1e-12
        mu(idx) = 1;
    else
        a = grid(idx);
        b = grid(idx+1);

        mu(idx)   = (b-val)/(b-a);
        mu(idx+1) = (val-a)/(b-a);
    end
end
