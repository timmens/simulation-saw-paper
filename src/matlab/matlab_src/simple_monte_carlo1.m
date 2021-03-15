% simple monte carlo study that uses DGP 2 

profile on;
tic;

S    = [1, 2, 3];          
N    = 300;%[25, 50, 100, 200];
T    = 2 .^ 7 + 1;%2 .^ [5, 6, 7] + 1;
nSim = 500;
%nDGP = 7;
nS   = size(S, 2);
nN   = 1;%size(N, 2);
nT   = 1;%size(T, 2);

option.nGrid     = 50; 
option.maxLambda = 100; 
option.minLambda = 0.0001;



%dgp     = repmat((1:nDGP)', [nSim, 1]);
s_est   = nan(nT * nS * nN, 1);
mdcj    = nan(nT * nS * nN, 1);
mse     = nan(nT * nS * nN, 1);

rng(123) % set random number generator seed
parpool(4) % declare (local) parallel cluster with 4 workers 

maxIter = nT * nS * nN;

%for i = 1:nDGP
    for t = 1:nT
        t_tmp = T(t);
        for s = 1:nS
            s_tmp = S(s);
            [beta, tau] = make_beta(t_tmp, s_tmp);
            for n = 1:nN
                disp("n = " + n + "; s = " + s + "; t = " + t)
                
                n_tmp = N(n);
                
                s_est_tmp = 0;
                mdcj_tmp  = 0;
                mse_tmp   = 0;
                parfor r = 1:nSim
                %for r = 1:nSim
                    [Y, X]    = DGP(t_tmp, n_tmp, beta, 2);
                    
                    [tau_est,alpha,~,~,~,~,~] = panelpls(Y, X, n_tmp, option, 1);
                    
                    beta_est  = alpha2beta(alpha, tau_est);
                    
                    % s_est
                    s_est_tmp = s_est_tmp + (length(tau_est) - 2);
                    
                    % MDCJ
                    mdcj_tmp  = mdcj_tmp + MDCJ(tau, tau_est, S); 
                    
                    % MSE
                    mse_tmp   = mse_tmp  + mean((beta - beta_est).^2);
                end 
                
                index        = (t - 1) * nS * nN + (s - 1) * nN + n;
                s_est(index) = s_est_tmp / nSim;
                mdcj(index)  = mdcj_tmp / nSim;
                mse(index)   = mse_tmp / nSim;
                
                % where are we: 
                
                disp((index / maxIter) * 100 + " % done");
            end
        end
    end
%end

[cN, cS, cT] = ndgrid(N, S, T);
T = cT(:); S = cS(:); N = cN(:);
result_table = table(T, S, N, s_est, mse, mdcj); 
disp(result_table);
toc; 



file_name = "simulation-" + regexprep(regexprep(datestr(datetime), ' ','-'), ':', '-')+ ".csv";
writetable(result_table,file_name)

profile viewer
