function L=designObserver(A,C)
%DESIGNOBSERVER Designs a discrete-time Luenberger observer using LMI approach for discrete-time systems.
%
%   L = DESIGNOBSERVER(A,B) computes a state observer gain K by solving
%   Linear Matrix Inequalities (LMIs) using YALMIP. The approach guarantees
%   closed-loop system stability for given matrices.
%
%   Designs a single (common) observer gain L that stabilizes all local
%   models.
%
%   Inputs:
%       A - Cell array containing state matrices {A1, A2, ..., Ats}
%       C - matrix containing output matrix C
%
%   Output:
%       L - State observer gain matrix
%
%   Example usage:
%       L = DesignObserver({A1,A2},C);
%
%   Requirements:
%       - YALMIP toolbox installed
%       - LMI solver compatible with YALMIP (e.g., lmilab,sedumi)
%
%   See also:
%       sdpsettings, solvesdp, checkset

%% LMI CONTROLLER

%% set options for given solver
options=sdpsettings;
options.solver='sedumi';
options.showprogress=1;
options.verbose=1;
warning('off','YALMIP:strict');

%% size of matrix A
n = size(A{1},1);
m = size(C,1);
ts=size(A,1);

%% clear yalmip state
yalmip('clear');

%% define variable matrices
P=sdpvar(n,n,'symmetric');
N=sdpvar(n,m,'full');

%% define inequalities LMI
LMIs = [];

%% use LMI variables
for i=1:ts
main_lmi = [-P,         A{i}'*P'-N'*C';
            P*A{i}-C*N,   -P];

LMIs = [LMIs, P >= 1e-6*eye(n)];
LMIs = [LMIs, main_lmi <= -1e-6*eye(2*n)];
end

%% solve it
res = optimize(LMIs, [], options);

%%
if res.problem == 0
    fprintf('=================================================================================\n');
    fprintf('Problem ...: %d\nInfo ......: %s\n', res.problem, res.info);
    if all(checkset(LMIs)>=0)
        fprintf('Everythings is all right :-)\n');
        fprintf('=================================================================================\n');
    else
        %disp(LMIs);
        %checkset(LMIs)
        fprintf('### Some LMIs are not feasible :-( ###\n');
        fprintf('=================================================================================\n');
    end
    rP = double(P);
    rN=double(N);
    L = rP\rN;
else
    res.problem
    L=ones(1,n);
end

% end