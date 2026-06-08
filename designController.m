function K=designController(A,B)
%DESIGNCONTROLLER Designs a state feedback discrete-time controller using LMI approach.
%
%   K = DesignController(A,B)
%
%   K = DESIGNCONTROLLER(A,B) computes a state feedback gain K by solving
%   Linear Matrix Inequalities (LMIs) using YALMIP. The approach guarantees
%   closed-loop system stability for given matrices.
%
%   Designs a single (common) state-feedback gain K that stabilizes all local
%   discrete-time models in {A{i},B{i}} (polytopic/TS vertex set).
%
%   Inputs:
%       A - Cell array containing state matrices {A1, A2, ..., Ats}
%       B - Cell array containing input matrices {B1, B2, ..., Bts}
%
%   Output:
%       K - State feedback controller gain matrix
%
%   Example usage:
%       K = DesignController({A1,A2},{B1,B2});
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

%% size of matrices
n = size(A{1},1);
r = size(B{1},2);
ts = size(A,1);

%% clear yalmip state
yalmip('clear');

%% define variable matrices
P=sdpvar(n,n,'symmetric');
N=sdpvar(r,n,'full');

%% define inequalities LMI
LMIs = [];
K=cell(1,ts);

%% use LMI variables
for i=1:ts
main_lmi = [-P,         P'*A{i}'-N'*B{i}';
            A{i}*P-B{i}*N,   -P];

LMIs= [LMIs, P >= 1e-6*eye(n)];
LMIs= [LMIs, main_lmi <= -1e-6*eye(2*n)];
end

%% solve it
res=optimize(LMIs,[], options);

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
    K = rN/rP;
else
    res.problem
    K=ones(1,n);
end

% end