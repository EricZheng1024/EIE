function [Sigma, exitflag] = update_CMAES(Sigma,X,lambda_def,flag_inj)
%Update the parameters in CMA-ES
% Remarks:
%   "h_sigma" is more likely to be true during the later stages of optimization.
%   "gen" in calculating "NoEffectAxis" is used for rotating indices.

    n = size(X,2);
    Sigma.gen = Sigma.gen + 1;

    % Initialize PSA
    if isempty(Sigma.PSA_buff)
        Sigma.PSA_buff.gamma_sigma_t = 0;
        Sigma.PSA_buff.gamma_c_t = 0;
        Sigma.PSA_buff.p_theta_t = zeros(n+(1+n)*n/2,1);
        Sigma.PSA_buff.gamma_theta_t = 0;
        Sigma.PSA_buff.lambda_t = lambda_def;
        Sigma.PSA_buff.sigma_star_t = [];
    end
    gamma_sigma_t = Sigma.PSA_buff.gamma_sigma_t;
    gamma_c_t = Sigma.PSA_buff.gamma_c_t;

    % Calculate the CMA parameters
    % Selection
    mu              = floor(Sigma.lambda/2);
    w               = log((Sigma.lambda+1)/2) - log(1:mu);
    w               = w./sum(w);  % the sum of the weights equals 1
    mu_eff          = (sum(w)^2) / sum(w.^2);
    c_m             = 1;
    % Injection
    cy              = sqrt(n) + 2*n/(n+2);  
    delta_max_sigma = 1;
    % Adaptation
    c_c             = (4+mu_eff/n) / (n+4+2*mu_eff/n);
    c_sigma         = (mu_eff+2) / (n+mu_eff+5);
    c_1             = 2 / ((n+1.3)^2+mu_eff);
    c_mu            = min(1-c_1,2*(mu_eff-2+1/mu_eff) / ((n+2)^2+mu_eff));
    d_sigma         = 1 + 2*max(0,sqrt((mu_eff-1)/(n+1))-1) + c_sigma;
    ENI             = sqrt(n)*(1-1/4/n+1/21/n^2);
    
    % Update the CMA model
    % Selection and recombination
    y               = (X(1:mu,:)-Sigma.m)/Sigma.sigma;  % row vector
    flag_inj        = flag_inj(1:mu);
    % y(flag_inj,:)   = min(1, cy./norm(Sigma.C^(-1/2)*y(flag_inj,:)')) * y(flag_inj,:);  % injection
    y(flag_inj,:)   = min(1, cy./vecnorm(Sigma.C^(-1/2)*y(flag_inj,:)',2,1))' .* y(flag_inj,:);  % injection
    y_w             = w*y;
    m_t             = Sigma.m;
    Sigma.m         = Sigma.m + c_m*Sigma.sigma*y_w;
    % Step-size control
    sigma_t         = Sigma.sigma;
    Sigma.p_sigma   = (1-c_sigma)*Sigma.p_sigma + sqrt(c_sigma*(2-c_sigma)*mu_eff)*Sigma.C^(-1/2)*y_w';
    gamma_sigma_t1  = (1-c_sigma)^2*gamma_sigma_t + c_sigma*(2-c_sigma);  % PSA

    % Sigma.sigma     = Sigma.sigma*exp(min(delta_max_sigma,c_sigma/d_sigma*(norm(Sigma.p_sigma)/ENI-1)));  % injection
    Sigma.sigma     = Sigma.sigma*exp(min(delta_max_sigma,c_sigma/d_sigma*(norm(Sigma.p_sigma)/ENI-sqrt(gamma_sigma_t1))));  % injection  PSA

    % Covariance matrix adaptation
    C_t             = Sigma.C;

    % h_sigma         = norm(Sigma.p_sigma)./sqrt(1-(1-c_sigma).^(2*Sigma.gen)) < (1.4+2/(n+1))*ENI;
    h_sigma         = norm(Sigma.p_sigma)./sqrt(1-(1-c_sigma).^(2*Sigma.gen)) < (1.4+2/(n+1))*ENI*sqrt(gamma_sigma_t1);  % PSA

    Sigma.p_c       = (1-c_c)*Sigma.p_c + h_sigma*sqrt(c_c*(2-c_c)*mu_eff)*y_w;
    w_circ          = w;  % all weights are required to be non-negative
    gamma_c_t1      = (1-c_c)^2*gamma_c_t + h_sigma*c_c*(2-c_c);  % PSA

    % Sigma.C         = (1+c_1*((1-h_sigma)*c_c*(2-c_c))-c_1-c_mu*sum(w,2))*Sigma.C + c_1*Sigma.p_c'*Sigma.p_c + c_mu*y'*diag(w_circ)*y;
    Sigma.C         = (1-c_1*gamma_c_t1-c_mu*sum(w,2))*Sigma.C + c_1*Sigma.p_c'*Sigma.p_c + c_mu*y'*diag(w_circ)*y;  % PSA

    Sigma.C         = triu(Sigma.C) + triu(Sigma.C,1)'; % Enforce symmetry

    % PSA
    [Sigma.lambda, Sigma.sigma, Sigma.PSA_buff] = PSA(Sigma.m, m_t, C_t, sigma_t, Sigma.C, Sigma.sigma, gamma_sigma_t1, gamma_c_t1, ...
        n, lambda_def, mu_eff, c_sigma, d_sigma, c_mu, c_c, c_1, ENI, w, X(1:mu,:), Sigma.PSA_buff);
    Sigma.PSA_buff.gamma_sigma_t = gamma_sigma_t1;
    Sigma.PSA_buff.gamma_c_t = gamma_c_t1;
    
    Sigma.C = triu(Sigma.C) + triu(Sigma.C,1)'; % Enforce symmetry
    
    % Reset the CMA model if possible
    exitflag = 0;
    [B,D] = eig(Sigma.C);  % The matrix B contains right eigenvectors as its columns, while D represents a diagonal matrix consisting of corresponding eigenvalues.
    diagD = diag(D);
    diagC = diag(Sigma.C);
    % Routine
    % NoEffectAxis  = all(Sigma.m==Sigma.m+0.1*Sigma.sigma*sqrt(diagD(mod(Sigma.gen,n)+1))*B(mod(Sigma.gen,n)+1,:));
    NoEffectAxis  = all(Sigma.m==Sigma.m+0.1*Sigma.sigma*sqrt(diagD).*B, 'all');
    % NoEffectCoord = all(Sigma.m==Sigma.m+0.2*Sigma.sigma*diagC');
    NoEffectCoord = all(Sigma.m==Sigma.m+0.2*Sigma.sigma*sqrt(diagC)');
    TolFun        = Sigma.gen_TolFun >= 10 + ceil(30*n/Sigma.lambda);
    % TolX          = all(Sigma.sigma*sqrt(diagC) < Sigma.sigma_0*1e-12) && all(Sigma.sigma*Sigma.p_c < Sigma.sigma_0*1e-12);  % default
    TolX          = all(Sigma.sigma*sqrt(diagC) < Sigma.sigma_0*1e-6) && all(Sigma.sigma*Sigma.p_c < Sigma.sigma_0*1e-6);
    % Exception
    % ConditionCov    = cond(Sigma.C) > 1e14;  % Condition number    If this termination condition is not included, the use of PSA results in "ecmnfish reporting a matrix singularity warning.
    ConditionCov    = false;
    TolXUp          = any(Sigma.sigma*sqrt(diagD) > 1e4*Sigma.sigma_0*sqrt(Sigma.diagD_0));
    ConditionNan    = isnan(Sigma.lambda) || any(isnan(Sigma.m)) || isnan(Sigma.sigma) || any(isnan(Sigma.C),'all') || any(isnan(Sigma.p_c)) || any(isnan(Sigma.p_sigma));

    if NoEffectAxis || NoEffectCoord || (TolFun && TolX)
        exitflag = -1;
    elseif ConditionCov || TolXUp || ConditionNan
        exitflag = 1;
    end
    if exitflag ~= 0
        Sigma = struct('s',[],'lambda',lambda_def,'m',[],'sigma',[],'sigma_0',[],'C',[],'diagD_0',[],'p_c',0,'p_sigma',0,'gen',0,'gen_TolFun',0,'PSA_buff',[]);
    end
end


function [lmabda_r_t1, sigma_t1, PSA_buff] = PSA(m_t1, m_t, C_t, sigma_t, C_t1, sigma_t1, gamma_sigma_t1, gamma_c_t1, ...
    n, lambda_def, mu_eff, c_sigma, d_sigma, c_mu, c_c, c_1, ENI, w, X, PSA_buff)
%Population size adaptation for CMA-ES
% ref: https://github.com/crucis/CPC881-PSA-CMA-ES

    p_theta_t = PSA_buff.p_theta_t;
    gamma_theta_t = PSA_buff.gamma_theta_t;
    lambda_t = PSA_buff.lambda_t;
    sigma_star_t = PSA_buff.sigma_star_t;

    alpha = 1.4;
    beta = 0.4;
    lambda_min = lambda_def;
    % lambda_max = inf;  % ref: "PSA-CMA-ES: ..."
    % lambda_max = 512*lambda_def;  % ref: "Benchmarking the PSA-CMA-ES on the BBOB noiseless testbed"
    lambda_max = 8*lambda_def;  % my configuration


    % Square root of Fisher information matrix    ref: https://github.com/crucis/CPC881-PSA-CMA-ES
    try
        sqrt_FIM = sqrtm(ecmnfish(X, C_t1));
    catch
        lmabda_r_t1 = nan;
        return
    end

    % Delta theta
    Delta_Sigma = sigma_t1^2 * C_t1 - sigma_t^2 * C_t;

    % tril_Delta_Sigma = tril(Delta_Sigma);  % The matrix is upper triangular in the paper, but the reproduced code employs a lower triangular matrix.
    % vech_Delta_Sigma = tril_Delta_Sigma(tril_Delta_Sigma~=0);  % If the matrix contains zeros, this may be trouble.

    vech_Delta_Sigma = Delta_Sigma( triu(true(size(Delta_Sigma))) );

    Delta_m = m_t1 - m_t;
    Delta_theta = [Delta_m, vech_Delta_Sigma']';

    % Expected value (denominator)
    first_bracket = 1 + 8 * gamma_sigma_t1 * (n - ENI^2)/ENI^2 * (c_sigma/d_sigma)^2;
    second_bracket = (n^2+n)*c_mu^2/mu_eff + (n^2+n)*c_c*(2-c_c)*c_1*c_mu*mu_eff*sum(w.^3) + c_1^2*(gamma_c_t1^2*n^2 + (1-2*gamma_c_t1+2*gamma_c_t1^2)*n);
    E_FIM_theta = (n * c_mu^2)/mu_eff + (2*n*(n-ENI^2)/ENI^2)*gamma_sigma_t1*(c_sigma/d_sigma)^2 + 0.5*first_bracket*second_bracket;


    % Update evolution path and its factor
    p_theta_t1 = (1-beta)*p_theta_t + sqrt(beta*(2-beta)) * (sqrt_FIM*Delta_theta)/sqrt(E_FIM_theta);
    gamma_theta_t1 = (1-beta)^2*gamma_theta_t + beta*(2-beta);

    % Update population size
    lambda_t1 = lambda_t * exp(beta*(gamma_theta_t1-norm(p_theta_t1)^2/alpha));
    lambda_t1 = min(max(lambda_t1,lambda_min),lambda_max);
    lmabda_r_t1 = round(lambda_t1);

    % Step-size correction
    % c_t1 = -sum(w'.*mean(X,2));  % This may be a bug in the reproduced code, which may result in sigma diverging or becoming negative.
    mu = floor(lmabda_r_t1/2);
    w = log((lmabda_r_t1+1)/2) - log(1:mu);
    w = w./sum(w);  % the sum of the weights equals 1
    c_t1 = -  w  *  norminv( ((1:mu)-0.375) / ( lmabda_r_t1-2*0.375+1) )';
    sigma_star_t1 = c_t1*n*mu_eff / (n-1+c_t1^2*mu_eff);
    if ~isempty(sigma_star_t)
        sigma_t1 = sigma_t1 * sigma_star_t1 / sigma_star_t;
    end

    % Update buff
    PSA_buff.p_theta_t = p_theta_t1;
    PSA_buff.gamma_theta_t = gamma_theta_t1;
    PSA_buff.lambda_t = lambda_t1;
    PSA_buff.sigma_star_t = sigma_star_t1;
end
