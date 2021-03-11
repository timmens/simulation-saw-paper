nSim    = 5000;
T       = 100;
N       = 100;
sigma_u = 1;

beta = zeros(T, 1);
beta((T/2 + 1):T) = 1;
beta = repmat(beta, [N, 1]);

option.nGrid     = 50; 
option.maxLambda = 100; 
option.minLambda = 0.0001;

rng(123)

s_est_tmp = 0;
parpool(4)
parfor r = 1:nSim
    % generate data
    xit = normrnd(0, 1, [T, N]);
    mui = mean(xit)';
    mui = repelem(mui, T);
    xit = reshape(xit, [T * N, 1]);
    uit = normrnd(0, 1, [T * N, 1]);
    
    yit = beta .* xit + mui + sigma_u * uit;
    
    % estimate
    [tau_est,~,~,~,~,~,~] = panelpls(yit, xit, N, option, 1);
                    
    %beta_est  = alpha2beta(alpha, tau_est);
                    
    s_est_tmp = s_est_tmp + length(tau_est); 
    
    % where am i 
    disp((r / nSim) * 100 + " % done");
end
s_est = s_est_tmp / nSim;
s_table = table(s_est);
writetable(s_table, "reconstruction" + regexprep(datestr(datetime), ' ', '-')+ ".csv")