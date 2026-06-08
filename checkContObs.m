function [co,ob] = checkContObs(Ac,Bc,C)
%CHECKCONTOBS Checks the controllability and observability of TS models.
%
%   [CO,OB] = CHECKCONTOBS(AC,BC,C) verifies the controllability
%   and observability of multiple Takagi-Sugeno (TS) models.
%
%   Inputs:
%       AC - Cell array containing state matrices for each TS model
%       BC - Cell array containing input matrices corresponding to AC
%       C  - Output matrix (assumed common for all models)
%
%   Outputs:
%       CO - 
%       OB - Logical vector indicating observability of each model
%
%   Example usage:Logical vector indicating controllability of each model
%       [co, ob] = checkContObs(Ac, Bc, C);
%
%   Note:
%       A warning will be displayed if any model is found not
%       controllable or observable.

nts = numel(Ac);
co = false(1, nts);
ob = false(1, nts);

for i = 1:nts
    co(i) = rank(ctrb(Ac{i}, Bc{i})) == size(Ac{i},1);
    ob(i) = rank(obsv(Ac{i}, C)) == size(Ac{i},1);
end

if all(co)
    disp('All TS models are controllable');
else
    warning('Some TS models are NOT controllable. Check "co" vector.');
end

if all(ob)
    disp('All TS models are observable');
else
    warning('Some TS models are NOT observable. Check "ob" vector.');
end
