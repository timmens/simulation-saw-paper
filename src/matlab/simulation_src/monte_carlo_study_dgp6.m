tic;
    
test = false;
unix = false;

N    = [30, 60, 120, 300];
T    = 2 .^ [5, 6, 7] + 1;
nSim = 500;
nN   = size(N, 2);
nT   = size(T, 2);

option.nGrid     = 20; 
option.maxLambda = 100; 
option.minLambda = 0.0001;

s_est_mean = nan(nT * nN, 1);
s_est_sd   = nan(nT * nN, 1);
mise_mean  = nan(nT * nN, 1);
mise_sd    = nan(nT * nN, 1);
s_0        = nan(nT * nN, 1);
taed_mean  = nan(nT * nN, 1);
taed_sd    = nan(nT * nN, 1);

rng_number = 123;
rng(rng_number) % set random number generator seed
parpool(20) % declare (local) parallel cluster with 4 workers 

nIter = nT * nN;

for t = 1:nT
    t_tmp = T(t);
    for n = 1:nN
        n_tmp = N(n);
        disp("dgp = " + 6 + "; n = " + n_tmp + "; t = " + t_tmp)

        s_est_mean_tmp = 0;
        s_est_sd_tmp = 0;
        mdcj_mean_tmp  = 0;
        mdcj_sd_tmp = 0;
        mise_mean_tmp  = 0;
        mise_sd_tmp = 0;
        hd_mean_tmp = 0;
        hd_sd_tmp = 0;
        s_0_tmp   = 0;
        taed_mean_tmp = 0;
        taed_sd_tmp = 0;

        parfor r = 1:nSim
            [Y, X, beta] = dgp6(t_tmp, n_tmp);

            [tau_est,alpha_est,~,~,~,~,~] = panelpls(Y, X, n_tmp, option, 1);

            beta_est = alpha2beta(alpha_est, tau_est);

            tau_est = tau_est(2:end-1); % they also report first and last time index
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

            % MISE mean
            mise_mean_tmp = mise_mean_tmp + mean((beta - beta_est).^2);
            % MISE sd
            mise_sd_tmp = mise_sd_tmp + mean((beta - beta_est).^2)^2;
            
            % Time-average euclidian distance
            gamma_true = beta_to_gamma(beta, []);
            gamma = beta_to_gamma(beta_est, []);
            taed_ = dist_euclidian_time_average(gamma, gamma_true);
            taed_mean_tmp = taed_mean_tmp + taed_;
            taed_sd_tmp = taed_sd_tmp + taed_^2;

            % Check if no jump location is found
            if isempty(tau_est)
                s_0_tmp = s_0_tmp + 1;
            end
        end

        index = (t - 1) * nN + n;
        s_est_mean(index) = s_est_mean_tmp  / nSim;
        s_est_sd(index)   = sqrt((s_est_sd_tmp - s_est_mean_tmp^2 / nSim) / (nSim - 1));
        mise_mean(index)  = mise_mean_tmp / nSim;
        mise_sd(index) = sqrt((mise_sd_tmp - mise_mean_tmp^2 / nSim) / (nSim - 1));
        s_0(index)   = s_0_tmp / nSim;
        taed_mean(index) = taed_mean_tmp / nSim;
        taed_sd(index) = sqrt((taed_sd_tmp - taed_mean_tmp^2 / nSim) / (nSim - 1));

        % where are we:
        disp((index / nIter) * 100 + " % done");
    end
end

[cN, cT] = ndgrid(N, T);
T = cT(:); N = cN(:);
result_table = table(T, N, s_est_mean, s_est_sd, mise_mean, mise_sd, s_0, taed_mean, taed_sd);

% getting bld path (matlab is such a *** language)
current_file_path = matlab.desktop.editor.getActiveFilename;
if unix
    splitter = "/src/";
    bld_suffix = "/bld/matlab/";
else
    splitter = "\src\";
    bld_suffix = "\bld\matlab\";
end

bld = split(current_file_path, splitter);
bld = string(bld);
bld = bld(1);
bld = bld + bld_suffix;

file_name = "simulation_dgp6";
if test
    file_name = file_name + "_" + regexprep(regexprep(datestr(datetime), ' ','-'), ':', '-')+ ".csv";
else
    file_name = file_name + ".csv";
end
writetable(result_table, bld + file_name)

% save additional information
ellapsed_time = toc / 60;
additional_info        = string(nan(4, 1));
additional_info(1, 1)  = "nsim = " + nSim;
additional_info(2, 1)  = "rng = " + rng_number;
additional_info(3, 1)  = "n.grid = " + option.nGrid;
additional_info(4, 1)  = "ellapsed time = " + ellapsed_time;

fid = fopen(bld + "additional-info_dgp6.txt", "w");
fprintf(fid, "%s\n", additional_info{:});
fclose(fid);
