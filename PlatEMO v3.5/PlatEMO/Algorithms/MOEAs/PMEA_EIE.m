classdef PMEA_EIE < ALGORITHM
% <multi/many> <real>
% PMEA
% tol --- 0.05 --- 
% is_EIE --- 1 --- 

%------------------------------- Reference --------------------------------
% "Solving Many-Objective Optimization Problems by a Pareto-Based
% Evolutionary Algorithm With Preprocessing and a Penalty Mechanism"
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [tol, is_EIE] = Algorithm.ParameterSet(0.05, 1);
            if length(tol) == 1
                tol = tol*ones(1,Problem.M);
            end

            %% Init MOEA
            % Generate the random population
            Population = Problem.Initialization();
            [zmin,I] = min(Population.objs,[],1);

            %% Init EIE
            % Parameters of CMA-ES
            if is_EIE
                state_CMAES = ones(1, Problem.M);  % ==1: activated; >=2: deactivated
            else
                state_CMAES = inf(1, Problem.M);
            end
            lambda_def = 4+floor(3*log(Problem.D));
            x_range = Problem.upper-Problem.lower;
            Sigma = struct('s',num2cell(1:Problem.M),'lambda',lambda_def,'m',num2cell(Population(I).decs,2)', ...
                'sigma',0.5*norm(diag(x_range)),'sigma_0',0.5*norm(diag(x_range)), ...  % 不能为0，会使得下次更新出现nan值
                'C',diag(x_range/norm(diag(x_range))),'diagD_0',(x_range/norm(diag(x_range))'), ...  % 考虑变量量纲不一致，并确保该对角阵二范数等于单位阵的二范数
                'p_c',0,'p_sigma',0,'gen',0,'gen_TolFun',0,'PSA_buff',[]);

            % Warm-start CMA-ES
            alpha = (tol*Problem.M) ./ (Problem.M-1-tol*(1-Problem.M));  % The meaning of alpha in the code differs slightly from that described in the paper; however, it is equivalent.
            WS_gamma = 0.1*ones(1,Problem.M);
            WS_alpha = 0.1*ones(1,Problem.M);
            objs_pop = Population.objs;
            objs_pop = normalize_2(objs_pop, zmin, max(Population.objs,[],1));
            objs_pop_trans = g_ews(objs_pop, objs_pop, alpha);
            for i = 1 : Problem.M
                [Sigma(i).m, Sigma(i).sigma, Sigma(i).C] = warm_start([Population.decs, objs_pop_trans(:,i)], WS_gamma(i), WS_alpha(i));
                Sigma(i).sigma_0 = Sigma(i).sigma;
                [~,D] = eig(Sigma(i).C);
                Sigma(i).diagD_0 = diag(D);
            end
   
            %% Optimization
            while Algorithm.NotTerminated(Population)
                % Generate offspring by EIE
                Offspring = [];
                Combine = cell(1, Problem.M);
                flag_inj = cell(1, Problem.M);
                for i = 1 : Problem.M
                    if state_CMAES(i) < 2
                        points_sampling = mvnrnd(Sigma(i).m,Sigma(i).sigma^2*Sigma(i).C,Sigma(i).lambda);
                        flag_inj{i} = any(points_sampling < Problem.lower, 2) | any(points_sampling > Problem.upper, 2);  % box constraints
                        Offspring_i = SOLUTION( min(max(points_sampling,Problem.lower),Problem.upper) );
                        Offspring = [Offspring Offspring_i];
                        Combine{i} = Offspring_i;
                    end
                end


                % Iteration of MOEA
                Offspring  = [Offspring, ...
                    OperatorDE(Population(randi(Problem.N,[1,Problem.N])),Population(randi(Problem.N,[1,Problem.N])),Population(randi(Problem.N,[1,Problem.N])),{0.9,0.5,1,50})];
                zmax = max(Population.objs,[],1);
                zmin = min([zmin;Offspring.objs],[],1);
                Population = EnvironmentalSelection([Population,Offspring],Problem.N, zmin);


                % Update EIE
                objs_off = Offspring.objs;
                objs_off = normalize_2(objs_off, zmin, zmax);
                objs_off_trans = g_ews(objs_off, objs_off, alpha);
                [~, I] = sort(objs_off_trans, 1);
                for i = 1 : Problem.M
                    % Fitness value stagnation detection
                    if state_CMAES(i) < 2
                        if ~exist('sol_fitness_best', 'var')
                            sol_fitness_best = Offspring(I(1,:));
                        else
                            objs_sfb = sol_fitness_best(i).obj;
                            objs_sfb = normalize_2(objs_sfb, zmin, zmax);
                            objs_sfb_trans_part = g_ews(objs_sfb(i),objs_sfb,alpha(i));
                            if objs_sfb_trans_part > objs_off_trans(I(1,i),i)
                                sol_fitness_best(i) = Offspring(I(1,i));
                                objs_sfb = sol_fitness_best(i).obj;
                                objs_sfb = normalize_2(objs_sfb, zmin, zmax);
                                objs_sfb_trans_part = g_ews(objs_sfb(i),objs_sfb,alpha(i));
                            end
                            objs_comb = Combine{i}.objs;
                            objs_comb = normalize_2(objs_comb, zmin, zmax);
                            diff_comb_sfb = g_ews(objs_comb(:,i), objs_comb, alpha(i)) - objs_sfb_trans_part;
                            if all(diff_comb_sfb < 1e-3)  % default 1e-12
                                Sigma(i).gen_TolFun = Sigma(i).gen_TolFun + 1;
                            else
                                Sigma(i).gen_TolFun = 0;
                            end
                        end
                    end
                    % Inject good solutions from the other new solutions
                    if state_CMAES(i) < 2 && Sigma(i).lambda == lambda_def  % 首先快速收敛，之后为了避免陷入局部最优，不允许injection
                        objs_comb = Combine{i}.objs;
                        [~,ia_objs_comb,ib_objs_off_top] = intersect(objs_comb,objs_off(I(1:Sigma(i).lambda,i),:),'stable','rows');
                        if length(ia_objs_comb) < Sigma(i).lambda  % 不是该CMA-ES产生的好解的逻辑索引存在1
                            index_objs_off_top = I(1:Sigma(i).lambda,i);
                            Combine{i} = [Combine{i}(ia_objs_comb) Offspring(index_objs_off_top(setdiff(1:Sigma(i).lambda,ib_objs_off_top,'stable')))];
                            flag_inj{i} = [flag_inj{i}(ia_objs_comb); true(Sigma(i).lambda-length(ib_objs_off_top),1)];
                        end
                    end
                end
                for i = 1 : Problem.M
                    if state_CMAES(i) < 2
                        % Sort the offsprings and injected solutions
                        objs_comb = Combine{i}.objs;
                        objs_comb = normalize_2(objs_comb, zmin, zmax);
                        objs_comb_trans_part = g_ews(objs_comb(:,i), objs_comb, alpha(i));
                        [~,rank] = sort(objs_comb_trans_part);
                        % Update CMA-ES
                        [Sigma(i), exitflag] = update_CMAES(Sigma(i),Combine{i}(rank).decs,lambda_def,flag_inj{i}(rank));
                        switch exitflag
                            case 0
                            case {-1, 1}
                                switch exitflag
                                    case -1
                                        state_CMAES(i) = state_CMAES(i) + 1;
                                    case 1
                                        % need restart
                                end
                                Sigma(i).s = i;
                                objs_pop = Population.objs;
                                objs_pop = normalize_2(objs_pop, zmin, zmax);
                                objs_pop_trans_part = g_ews(objs_pop(:,i), objs_pop, alpha(i));
                                % Restart
                                % [~, I] = min(objs_pop_trans_part);
                                % Sigma(i).m = Population(I).dec;
                                % x_range = Problem.upper-Problem.lower;
                                % Sigma(i).sigma = 0.5*norm(diag(x_range));
                                % Sigma(i).C = diag(x_range/norm(diag(x_range)));
                                % Warm restart
                                [Sigma(i).m, Sigma(i).sigma, Sigma(i).C] = warm_start([Population.decs, objs_pop_trans_part], WS_gamma(i), WS_alpha(i));
                                Sigma(i).sigma_0 = Sigma(i).sigma;
                                [~,D] = eig(Sigma(i).C);
                                Sigma(i).diagD_0 = diag(D);
                        end
                    end
                end
            end
        end
    end
end


%%
function Population = EnvironmentalSelection(Population,N,Zmin)
% The environmental selection of PMEA

    [FrontNo,MaxFNo] = NDSort(Population.objs,N);
    St = find(FrontNo <= MaxFNo);
    
    % Selection
    if length(St) == N
        Population = Population(St);
    else
        Population = LastSelection(Population(St), N, Zmin);
    end
end


function Population = LastSelection(Population, K, Zmin)

    PopObj = Population.objs;
    Index = true(1,size(PopObj,1));
    Fitness = sum(PopObj,2);
    Fit = sort(Fitness);
    if ceil(length(Fit)*3/4) < K    %如果总数的0.75比N小，也就是删除0.25的解之后数量到不了N
        q = Fit(K)+1.5*(Fit(K) - Fit(size(PopObj,1)-K));
    else
        q = Fit(ceil(length(Fit)*3/4))+1.5*(Fit(ceil(length(Fit)*3/4)) - Fit(ceil(length(Fit)/4)));
    end
    %q是上四分位点的界
    c = find(Fitness(Index) > q);%找到大于上四分位点的解
    if ~isempty(c)
        Index(c) = false;%如果存在大于四分位点的数，将这些数删掉
    end
    %只有四分之三的解称为有效解
    D = sum(max(PopObj(Index,:),[],1)-Zmin);        %剩下的最大的解（有效解）减去参考点之后
    Is = sum(max(PopObj,[],1)-Zmin) / D;    %从所有的种群中找到最大值与之前的0.75中的最大值进行比较
    
    index = find(Index==false);
    %找到多出来的那0.25
    if Is > 1.1  %大于1.1，说明这部分多出来的解是DRS
        PopObj = PopObj(Index,:);
        Population = Population(Index);
    end
    Zmax = max(PopObj,[],1);%最大的参考点
    [N,M]  = size(PopObj);%OBJ的规模
    PopObj = (PopObj-repmat(Zmin,N,1))./repmat(Zmax-Zmin,N,1);%归一化整个种群
    
    fitness = sum(PopObj,2);%将每个解的三个目标加在一起？
    angle = pdist2(PopObj,PopObj,'cosine');%
    angle(logical(eye(length(angle)))) = inf;
    %将自己和自己的距离定义为无穷
    
    a = max(min(angle,[],2)); %找到所有列中最小的数中最大的那一个
    pn = N-K;  %N为现在种群的实际大小，K是想要维持的种群的大小，pn是需要惩罚的数量
    Choose                          = false(1,N);
    zChoose                         = false(1,N);
    punish                          = false(1,N);
    if ~any(Choose)  %如果choose中都是0，这是1,出现一个1，就是0，any，哪怕其中有一个零的，any，都会给个1
        %         Select the extreme solutions first
        [~,extreme]        = min(pdist2(PopObj,eye(M),'cosine'),[],1);%eye，生成单位矩阵，extreme是M个极端解的索引
        Choose(extreme) = true;  % 将这三个解加入Choose中
        zChoose(extreme) = true;
        for i = 1 : M
            Have                             = find(~Choose); %Have 是除了极端解之外的索引值
            [~, pos]                   = find(angle(extreme(i),Have) < a);  %找极端解和非极端解的小于a的,小于a表示拥挤，返回位置
            Z = find(~punish(Have(pos(1:end))));  % 找到这几个解在整个种群的中的位置
            if sum(punish) + length(Z) <= pn   %惩罚集数量不足
                punish(Have(pos(1:end)))                  = true;  %将这些拥挤的解加入惩罚集
                zChoose(Have(pos(1:end)))                 = true;  %将这些拥挤解加入zchoose
            else
                Select = find(Choose);   %入股惩罚集数量超了
                punish(Have(pos(1:end))) = true;
                Punished = find(punish);
                [neib, ~] = min(angle(Punished,Select),[],2);
                [~,I] = sort(neib(1:end));
                zChoose(1:end) = false;
                punish(1:end) = false;
                punish(Punished(I(1:pn))) = true;
                zChoose([Select,Punished(I(1:pn))]) = true;   %将这些拥挤解加入ZChoose
            end
        end
    end
    while sum(Choose) < K    %当选中的解的数量小于K个预定解时
        Remain                           = find(~zChoose);  %是没有进入收敛集和拥挤集的
        [~, index]                       = min(fitness(Remain));  %寻找收敛性最好的第二个解
        Have                             = find(~Choose);   %have中是惩罚集中的解
        Choose(Remain(index))            = true;   %再次补充收敛性好的解
        zChoose(Remain(index))           = true;   %候补解状态置一
        [~, pos]                         = find(angle(Remain(index),Have) < a);%再找一下附近的拥挤解
        Z = find(~punish(Have(pos(1:end))));  % 需要再次惩罚的解
        if sum(punish) + length(Z) <= pn     %如果这次选出来的惩罚解满足数量要求
    
            punish(Have(pos(1:end)))                  = true;
            zChoose(Have(pos(1:end)))                 = true;
        else
            Select = find(Choose);   %
            punish(Have(pos(1:end))) = true;
            Punished = find(punish);
            [neib, ~] = min(angle(Punished,Select),[],2);
            [~,I] = sort(neib(1:end));
            zChoose(1:end) = false;
            punish(1:end) = false;
    
            punish(Punished(I(1:pn))) = true;
            zChoose([Select,Punished(I(1:pn))]) = true;
        end
    end
    Population = Population(Choose);
end

