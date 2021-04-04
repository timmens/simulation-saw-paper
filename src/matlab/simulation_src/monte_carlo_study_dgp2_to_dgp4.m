tic;

S    = [1, 2, 3];          
N    = [30, 60, 120, 300];
T    = 2 .^ [5, 6, 7] + 1;
dgp  = [2, 3, 4];

nSim = 500;
nDGP = size(dgp, 2);
nS   = size(S, 2);
nN   = size(N, 2);
nT   = size(T, 2);

option.nGrid     = 25;
option.maxLambda = 100;
option.minLambda = 0.0001;

s_est_mean = nan(nT * nS * nN * nDGP, 1);
s_est_sd   = nan(nT * nS * nN * nDGP, 1);
mdcj_mean  = nan(nT * nS * nN * nDGP, 1);
mise_mean  = nan(nT * nS * nN * nDGP, 1);
hd_mean    = nan(nT * nS * nN * nDGP, 1);
mdcj_sd    = nan(nT * nS * nN * nDGP, 1);
mise_sd    = nan(nT * nS * nN * nDGP, 1);
hd_sd      = nan(nT * nS * nN * nDGP, 1);
s_0        = nan(nT * nS * nN * nDGP, 1);


rng_number = 321;
rng(rng_number)
parpool(4) % declare (local) parallel cluster with 4 workers 

nIter = nT * nS * nN * nDGP;

for i = dgp
    for t = 1:nT
        t_tmp = T(t);
        for s = 1:nS
            s_tmp = S(s);
            [beta, tau] = make_beta(t_tmp, s_tmp);
            for n = 1:nN
                n_tmp = N(n);
                disp("DGP = " + i + "; n = " + n_tmp + "; s = " + s_tmp + "; t = " + t_tmp)
                
                s_est_mean_tmp = 0;
                s_est_sd_tmp = 0;
                mdcj_mean_tmp  = 0;
                mdcj_sd_tmp = 0;
                mise_mean_tmp  = 0;
                mise_sd_tmp = 0;
                hd_mean_tmp = 0;
                hd_sd_tmp = 0;
                s_0_tmp   = 0;

                parfor r = 1:nSim
                    if i == 2
                        [Y, X, Z] = dgp2(t_tmp, n_tmp, beta);
                        [tau_est, alpha_est] = panelpgmm(Y, X, Z, n_tmp, option);                        
                    else
                        [Y, X] = DGP(t_tmp, n_tmp, beta, i);
                        [tau_est,alpha_est] = panelpls(Y, X, n_tmp, option, 1);
                    end    
                    
                    beta_est  = alpha2beta(alpha_est, tau_est);
		    
                    tau_est   = tau_est(2:end-1); % they also report first and last time index
                    if ~isempty(tau_est)
                        tau_est = tau_est - 1;
                        % they report tau_i as starting of next regime we as end of current regime
                    end 
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%% CODE EXCEPTIONS %%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    % s_est mean
                    s_est_mean_tmp = s_est_mean_tmp + length(tau_est);
                    % s_est sd
                    s_est_sd_tmp = s_est_sd_tmp + length(tau_est)^2;
                    
                    % MDCJ mean
                    mdcj_mean_tmp = mdcj_mean_tmp + MDCJ(tau, tau_est, S);
                    % MDCJ sd
                    mdcj_sd_tmp = mdcj_sd_tmp + MDCJ(tau, tau_est, S)^2;
                    
                    % MISE mean
                    mise_mean_tmp = mise_mean_tmp + mean((beta - beta_est).^2);
                    % MISE sd
                    mise_sd_tmp = mise_sd_tmp + mean((beta - beta_est).^2)^2;
                    
                    % Hausdorff mean
                    hd_mean_tmp  = hd_mean_tmp + dist_hausdorff(tau, tau_est);
                    % Hausdorff sd
                    hd_sd_tmp  = hd_sd_tmp + dist_hausdorff(tau, tau_est)^2;

                    % Check if no jump location is found
                    if isempty(tau_est)
                        s_0_tmp = s_0_tmp + 1;
                    end

                end
                
                index = (i - 2) * nT * nS * nN + (t - 1) * nS * nN + (s - 1) * nN + n;

                s_est_mean(index) = s_est_mean_tmp  / nSim;
                s_est_sd(index)   = sqrt((s_est_sd_tmp - s_est_mean_tmp^2 / nSim) / (nSim - 1));
                mdcj_mean(index)  = mdcj_mean_tmp  / nSim;
                mdcj_sd(index) = sqrt((mdcj_sd_tmp - mdcj_mean_tmp^2 / nSim) / (nSim - 1));
                mise_mean(index)  = mise_mean_tmp / nSim;
                mise_sd(index) = sqrt((mise_sd_tmp - mise_mean_tmp^2 / nSim) / (nSim - 1));
                s_0(index)   = s_0_tmp / nSim;
                hd_mean(index)    = hd_mean_tmp / nSim;
                hd_sd(index) = sqrt((hd_sd_tmp - hd_mean_tmp^2 / nSim) / (nSim - 1));
                
                % where are we: 
                disp((index / nIter) * 100 + " % done");
            end
        end
    end
end

[cN, cS, cT, cDGP] = ndgrid(N, S, T, dgp);
T = cT(:); S = cS(:); N = cN(:); dgp = cDGP(:);
result_table = table(dgp, T, S, N, s_est_mean, s_est_sd, mise_mean, mise_sd, mdcj_mean, mdcj_sd, hd_mean, hd_sd, s_0); 

% getting bld path (matlab is such a *** language)
current_file_path = matlab.desktop.editor.getActiveFilename;
bld = split(current_file_path, "/src/");
bld = string(bld);
bld = bld(1);
bld = bld + "/bld/matlab/";

file_name = "simulation_dgp2-to-dgp4-" + regexprep(regexprep(datestr(datetime), ' ','-'), ':', '-')+ ".csv";
writetable(result_table, bld + file_name)

                s_est_mean(index) = s_est_mean_tmp  / nSim;
                s_est_sd(index)   = sqrt((s_est_sd_tmp - s_est_mean_tmp^2 / nSim) / (nSim - 1));
                mdcj_mean(index)  = mdcj_mean_tmp  / nSim;
                mdcj_sd(index) = sqrt((mdcj_sd_tmp - mdcj_mean_tmp^2 / nSim) / (nSim - 1));
                mise_mean(index)  = mise_mean_tmp / nSim;
                mise_sd(index) = sqrt((mise_sd_tmp - mise_mean_tmp^2 / nSim) / (nSim - 1));
                s_0(index)   = s_0_tmp / nSim;
                hd_mean(index)    = hd_mean_tmp / nSim;
                hd_sd(index) = sqrt((hd_sd_tmp - hd_mean_tmp^2 / nSim) / (nSim - 1));llapsed_time = toc; 
additional_info        = string(nan(4, 1));
additional_info(1, 1)  = "nsim = " + nSim;
additional_info(2, 1)  = "rng = " + rng_number;
additional_info(3, 1)  = "n.grid = " + option.nGrid;
additional_info(4, 1)  = "ellapsed time = " + ellapsed_time;

fid = fopen(bld + "additional-info_dgp2-to-dgp4.txt", "w");
fprintf(fid, "%s\n", additional_info{:});
fclose(fid);
