function objs_n = normalize_2(objs, lb, ub)
    objs_n = (objs - lb) ./ (ub - lb);
end