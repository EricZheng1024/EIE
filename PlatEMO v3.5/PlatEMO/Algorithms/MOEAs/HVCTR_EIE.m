classdef HVCTR_EIE < ALGORITHM
% <multi> <real>
% HVCTR
% r_small --- 1 --- reference point 1
% r_large --- 100 --- reference point 2
% tol --- 0.05 --- 
% is_EIE --- 1 --- 

%------------------------------- Reference --------------------------------
% "Hpervolume-based cooperative coevolution with two reference points for
% multi-objective optimization"
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [r_small, r_large, tol, is_EIE] = Algorithm.ParameterSet(1, 100, 0.05, 1);
            if length(tol) == 1
                tol = tol*ones(1,Problem.M);
            end

            %% Init MOEA
            % Generate random population
            Population1 = Problem.Initialization();
            Population2 = Problem.Initialization();

            N = Problem.N;
            N1 = floor(N /2); 
            N2 = N-N1;

            % DSS selection (without non-dominated sorting)
            [~, ind_sel1] = DSS(Population1.objs,N1);
            [~, ind_sel2] = DSS(Population2.objs,N2);

            Population1 = Population1(ind_sel1);     
            Population2 = Population2(ind_sel2);       
            Population = [Population1,Population2];

            % Initialize the ideal point
            [zmin, I] = min(Population.objs,[],1);

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
                % Update populaiton
                if ~isempty(Offspring)
                    zmin = min([zmin; Offspring.objs]);  % update the ideal point
                    Offspring = Offspring(randperm(length(Offspring)));
                    for i = 1 : length(Offspring)
                        Pt = [Population,Offspring(i)];
                        zmax = max(Pt.objs, [], 1);
                        [Population1] = EnvSel(Population1,Offspring(i),Population2,r_large,zmin,zmax);
                        [Population2] = EnvSel(Population2,Offspring(i),Population1,r_small,zmin,zmax);
                        Population = [Population1, Population2];
                    end
                end

                % Generate offspring by DE and update population
                for i = 1 : Problem.N
                    drawnow();

                    % Generate offspring1 from subpopulation1
                    Offspring1 = OperatorDE(Population1(randperm(end,1)),Population1(randperm(end,1)),Population1(randperm(end,1)),{0.9,0.5,1,50});
                    Offspring2 = OperatorDE(Population(randperm(end,1)),Population(randperm(end,1)),Population(randperm(end,1)),{0.9,0.5,1,50});
                    Offspring = [Offspring, Offspring1, Offspring2];

                    % Estimate Zmin and Zmax
                    AllPopObj = [Population1.objs; Offspring1.obj;Population2.objs;Offspring2.obj];
                    % Zmin      = min(AllPopObj,[],1);
                    zmin      = min([zmin; Offspring1.obj; Offspring2.obj], [], 1);
                    zmax      = max(AllPopObj,[],1);

                    % Update each subpopulation
                    [Population1] = EnvSel(Population1,Offspring(end-1:end),Population2,r_large,zmin,zmax);
                    [Population2] = EnvSel(Population2,Offspring(end-1:end),Population1,r_small,zmin,zmax);

                    % Merge two populations
                    Population = [Population1, Population2];
                end


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
function [Population] = EnvSel(Population,Offspring,Population_aux,r,Zmin,Zmax)
    for i = 1 : length(Offspring)
        Population = [Population, Offspring(i)];
        N     = length(Population);
        % Identify the solutions in the last front
        AllPopObj = [Population.objs;Population_aux.objs];
        [FrontNo,MaxFNo] = NDSort(AllPopObj,inf);
        while true
            LastFront = find(FrontNo==MaxFNo);
            % Include at least one solution from Population
            if sum(LastFront <= N) > 0
                break;
            end
            MaxFNo = MaxFNo - 1;
        end
    
        % Normalize the solutions in the last front
        AllPopObj_last    = AllPopObj(LastFront,:);
        AllPopObj_last_normalized = (AllPopObj_last - Zmin)./(Zmax-Zmin);
     
        [N_last,~]     = size(AllPopObj_last);
        
        % Calculate the contribution of hypervolume of each solution
        deltaS = CalHVC(AllPopObj_last_normalized,r,N_last);
    
        % Set HVC of aux Population to inf
        deltaS(LastFront > N) = inf;
        
        % Delete the worst solution from the last front
        [~,worst] = min(deltaS,[],'includenan');
        Population(LastFront(worst)) = [];      
    end
end

function HVC = CalHVC(data,r,PopNum)
 
    ref = ones(1,size(data,2)).*r;    
%     HVC = zeros(1,PopNum);
%     for i=1:PopNum
%         data1 = data;
%         s = data1(i,:);
%         data1(i,:)=[];
%         data1 = max(s,data1);        
%         HVC(1,i) = prod(ref-s)-stk_dominatedhv(data1,ref); 
%     end
    HV = Hypervolume_MEX(data,ref);
    
    if PopNum > 2
        hv = zeros(1,PopNum);
        hv(1) = Hypervolume_MEX(data(2:end,:),ref);
        hv(PopNum) = Hypervolume_MEX(data(1:end-1,:),ref);
        for i = 2 : PopNum-1
            Pop = [data(1:i-1,:);data(i+1:end,:)];
            hv(i) = Hypervolume_MEX(Pop,ref);
        end
        HVC = HV - hv;
    elseif PopNum < 2
        HVC = HV;
    else
        hv = zeros(1,PopNum);
        hv(1) = Hypervolume_MEX(data(end,:),ref);
        hv(PopNum) = Hypervolume_MEX(data(1,:),ref);
        HVC = HV - hv;
    end   
end


function [B, ind_sel] = DSS(a,sel_num)

if size(a,1)<=sel_num
    B=a; 
    ind_sel = 1: size(a,1);
else  
    %normalization 
    fmax=max(a);
    fmin=min(a);
    anorm=(a-fmin)./(fmax-fmin);
    
    [sz_a, M]=size(a);  
    [~, ind_extreme] = max(a(:,1));
    first_extreme=a(ind_extreme,:); %get the first extreme solution
    
    %initialize B as empty set 
    B=[]; 

    %initialize v and d 
    v=zeros(sz_a, 1);     %zeroes represent false 
    d=ones(sz_a,1)*inf;   %initiliaze distance to a very large number, i.e., inf

    %select first extreme solutions 
    B=[B;first_extreme]; 
    v(ind_extreme)=1; 

    for l=1:sz_a
        if v(l)==0
            d(l)=min(sqrt(sum((anorm(l,:)-anorm(ind_extreme,:)).^2)),d(l));
        end
    end

    sz_B=size(B,1); 
    ind_sel=[];
    ind_sel=[ind_sel;ind_extreme]; 
    %select an isolated solution to B 
    while sz_B < sel_num
        z=find(v==0); 
        max_distance=max(d(z)); 
        j=find(d==max_distance);
        j=j(1);
        B=[B;a(j,:)];
        v(j)=1; 
        d(j)=inf;

        for l=1:sz_a
            if v(l)==0
                d(l)=min(sqrt(sum((anorm(l,:)-anorm(j,:)).^2)),d(l));
            end
        end
        sz_B=size(B,1); 
        ind_sel=[ind_sel;j];
    end   
end
end
