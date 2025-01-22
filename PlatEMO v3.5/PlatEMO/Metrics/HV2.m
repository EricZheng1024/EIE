function score = HV2(Population,optimum)
% <max>
% Hypervolume
% 基于PlatEMO的HV.m修改了fmin、归一化、参考点
% 2022.8.29：修改解越界丢弃的条件，原为PopObj(any(PopObj>fmax,2),:)
% 2022.10.23：any(pop_min<fmin,2) || any(pop_min<0,2) || any(fmin<0,2)改为any(pop_min<fmin,2)
% 有些测试问题产生PF的方式是采用UniformPoint，边界的采样不会到0而是1e-6，因此HV2会返回nan

%------------------------------- Reference --------------------------------
% E. Zitzler and L. Thiele, Multiobjective evolutionary algorithms: A
% comparative case study and the strength Pareto approach, IEEE
% Transactions on Evolutionary Computation, 1999, 3(4): 257-271.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
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
        fmin = min(optimum,[],1);
        pop_min = min(PopObj,[],1);
        if any(pop_min<fmin,2)
            score = nan;
            return
        end
        fmax = max(optimum,[],1);
        PopObj(any(PopObj>(fmax+0.1*(fmax-fmin)),2),:) = [];
        [N,M] = size(PopObj);
        PopObj = (PopObj-repmat(fmin,N,1))./repmat((fmax-fmin),N,1);
        RefPoint = 1.1*ones(1,M);
        if isempty(PopObj)
            score = 0;
        elseif M < 4
            % Calculate the exact HV value
            pl = sortrows(PopObj);
            S  = {1,pl};
            for k = 1 : M-1
                S_ = {};
                for i = 1 : size(S,1)
                    Stemp = Slice(cell2mat(S(i,2)),k,RefPoint);
                    for j = 1 : size(Stemp,1)
                        temp(1) = {cell2mat(Stemp(j,1))*cell2mat(S(i,1))};
                        temp(2) = Stemp(j,2);
                        S_      = Add(temp,S_);
                    end
                end
                S = S_;
            end
            score = 0;
            for i = 1 : size(S,1)
                p     = Head(cell2mat(S(i,2)));
                score = score + cell2mat(S(i,1))*abs(p(M)-RefPoint(M));
            end
        else
            % Estimate the HV value by Monte Carlo estimation
            SampleNum = 1e6;
            MaxValue  = RefPoint;
            MinValue  = min(PopObj,[],1);
            Samples   = unifrnd(repmat(MinValue,SampleNum,1),repmat(MaxValue,SampleNum,1));
            for i = 1 : size(PopObj,1)
                drawnow();
                domi = true(size(Samples,1),1);
                m    = 1;
                while m <= M && any(domi)
                    domi = domi & PopObj(i,m) <= Samples(:,m);
                    m    = m + 1;
                end
                Samples(domi,:) = [];
            end
            score = prod(MaxValue-MinValue)*(1-size(Samples,1)/SampleNum);
        end
    end
end

function S = Slice(pl,k,RefPoint)
    p   = Head(pl);
    pl = Tail(pl);
    ql = [];
    S  = {};
    while ~isempty(pl)
        ql  = Insert(p,k+1,ql);
        p_  = Head(pl);
        cell_(1,1) = {abs(p(k)-p_(k))};
        cell_(1,2) = {ql};
        S   = Add(cell_,S);
        p   = p_;
        pl  = Tail(pl);
    end
    ql = Insert(p,k+1,ql);
    cell_(1,1) = {abs(p(k)-RefPoint(k))};
    cell_(1,2) = {ql};
    S  = Add(cell_,S);
end

function ql = Insert(p,k,pl)
    flag1 = 0;
    flag2 = 0;
    ql    = [];
    hp    = Head(pl);
    while ~isempty(pl) && hp(k) < p(k)
        ql = [ql;hp];
        pl = Tail(pl);
        hp = Head(pl);
    end
    ql = [ql;p];
    m  = length(p);
    while ~isempty(pl)
        q = Head(pl);
        for i = k : m
            if p(i) < q(i)
                flag1 = 1;
            else
                if p(i) > q(i)
                    flag2 = 1;
                end
            end
        end
        if ~(flag1 == 1 && flag2 == 0)
            ql = [ql;Head(pl)];
        end
        pl = Tail(pl);
    end  
end

function p = Head(pl)
    if isempty(pl)
        p = [];
    else
        p = pl(1,:);
    end
end

function ql = Tail(pl)
    if size(pl,1) < 2
        ql = [];
    else
        ql = pl(2:end,:);
    end
end

function S_ = Add(cell_,S)
    n = size(S,1);
    m = 0;
    for k = 1 : n
        if isequal(cell_(1,2),S(k,2))
            S(k,1) = {cell2mat(S(k,1))+cell2mat(cell_(1,1))};
            m = 1;
            break;
        end
    end
    if m == 0
        S(n+1,:) = cell_(1,:);
    end
    S_ = S;     
end
