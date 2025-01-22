function score = DistIdeal(Population,optimum,p)
% <min>
% Lp-distance of estimated ideal point to ideal point

%------------------------------- Reference --------------------------------
% 
%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    PopObj = Population.best.objs;
    if size(PopObj,2) ~= size(optimum,2)
        score = nan;
    else
        switch nargin
            case 2
                p = 2;
            case 3
                % do noting
            otherwise
                error('The number of input parameters is wrong.')
        end
        pop_min = min(PopObj, [], 1);
        opt_min = min(optimum, [], 1);
        opt_max = max(optimum, [], 1);
        score = norm((pop_min-opt_min)./(opt_max-opt_min), p);
    end
end