function res = g_ews(objs_part, objs, alpha)
%Calculate the EWS-based subproblem function value

    if size(alpha,2) ~= size(objs_part,2)
        error('Incorrect input parameters for cal_fitness().')
    end
    res = (1-alpha).*objs_part + alpha.*mean(objs,2);
end