% monte carlo study for DGP 1
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
s_est_sd = nan(nT * nN, 1);
mise1_mean = nan(nT * nN, 1);
mise2_mean = nan(nT * nN, 1);
mise1_sd = nan(nT * nN, 1);
mise2_sd = nan(nT * nN, 1);
hd1_mean = nan(nT * nN, 1);
hd2_mean = nan(nT * nN, 1);
hd1_sd = nan(nT * nN, 1);
hd2_sd = nan(nT * nN, 1);
taed_mean = nan(nT * nN, 1);
taed_sd = nan(nT * nN, 1);

rng_number = 321;
rng(rng_number) % set random number generator seed
parpool(20) % declare (local) parallel cluster with 4 workers 

nIter = nT * nN;

for t = 1:nT
    t_tmp = T(t);
    for n = 1:nN

        disp(strcat('dgp1: n = ', num2str(N(n)), '; t = ', num2str(T(t))));
        n_tmp = N(n);
        
        results = nan([10, nSim]);
        s_est_mean_tmp = 0;
        s_est_sd_tmp = 0;
        mise_mean_tmp1 = 0;
        mise_mean_tmp2 = 0;
        mise_sd_tmp1  = 0;
        mise_sd_tmp2  = 0;
        hd1_mean_tmp = 0;
        hd2_mean_tmp = 0;
        hd_sd_tmp1   = 0;
        hd_sd_tmp2   = 0;
        taed_mean_tmp = 0;
        taed_sd_tmp = 0;

        parfor r = 1:nSim
            [Y, X1, X2, tau1, tau2, beta1, beta2] = dgp1(t_tmp, n_tmp);

            [tau_est,alpha,~,~,~,~,~] = panelpls(Y, [X1, X2], n_tmp, option, 1);  % 1 for mex
            beta_est  = alpha2beta(alpha, tau_est);
            
            tau_est   = tau_est(2:end-1); % they also report first and last time index
            if ~isempty(tau_est)
                tau_est = tau_est - 1;
            end

            s_est_mean_tmp = s_est_mean_tmp + length(tau_est);
            s_est_sd_tmp = s_est_sd_tmp + length(tau_est)^2;

            % mise_mean
            mise_mean_tmp1 = mise_mean_tmp1 + mean((beta1 - beta_est(:, 1)).^2);
            mise_mean_tmp2 = mise_mean_tmp2 + mean((beta2 - beta_est(:, 2)).^2);
            % mise_sd
            mise_sd_tmp1 = mise_sd_tmp1 + mean((beta1 - beta_est(:, 1)).^2)^2;
            mise_sd_tmp2 = mise_sd_tmp2 + mean((beta2 - beta_est(:, 2)).^2)^2;
                                
            % Hausdorff
            hd1_mean_tmp = hd1_mean_tmp + dist_hausdorff(tau1, tau_est) / t_tmp;
            hd2_mean_tmp = hd2_mean_tmp + dist_hausdorff(tau2, tau_est) / t_tmp;
            % Hausdorff sd
            hd_sd_tmp1 = hd_sd_tmp1 + (dist_hausdorff(tau1, tau_est) / t_tmp)^2;
            hd_sd_tmp2 = hd_sd_tmp2 + (dist_hausdorff(tau2, tau_est) / t_tmp)^2;
            
            % Time-average of euclidian distance
            gamma_true = beta_to_gamma([beta1, beta2], []);
            gamma = beta_to_gamma(beta_est, []);
            taed_ = dist_euclidian_time_average(gamma, gamma_true);
            taed_mean_tmp = taed_mean_tmp + taed_;
            taed_sd_tmp = taed_sd_tmp + taed_^2;
            
        end 

        index         = (t - 1) * nN + n;
        s_est_mean(index)  = s_est_mean_tmp / nSim;
        s_est_sd(index) = sqrt((s_est_sd_tmp - s_est_mean_tmp^2 / nSim) / (nSim - 1));
        mise1_mean(index) = mise_mean_tmp1 / nSim;
        mise2_mean(index) = mise_mean_tmp2 / nSim;
        mise1_sd(index) = sqrt((mise_sd_tmp1 - mise_mean_tmp1^2 / nSim) / (nSim - 1));
        mise2_sd(index) = sqrt((mise_sd_tmp2 - mise_mean_tmp2^2 / nSim) / (nSim - 1));
        hd1_mean(index) = hd1_mean_tmp / nSim;
        hd2_mean(index) = hd2_mean_tmp / nSim;
        hd1_sd(index) = sqrt((hd_sd_tmp1 - hd1_mean_tmp^2 / nSim) / (nSim - 1));
        hd2_sd(index) = sqrt((hd_sd_tmp2 - hd2_mean_tmp^2 / nSim) / (nSim - 1));
        taed_mean(index) = taed_mean_tmp / nSim;
        taed_sd(index) = sqrt((taed_sd_tmp - taed_mean_tmp^2 / nSim) / (nSim - 1));

        % where are we: 
        disp(strcat(num2str((index / nIter) * 100), ' % done'));
    end
end


[cN, cT] = ndgrid(N, T);
T = cT(:); N = cN(:);
result_table = table(T, N, s_est_mean, s_est_sd, mise1_mean, mise1_sd, mise2_mean, mise2_sd, hd1_mean, hd1_sd, hd2_mean, hd2_sd, taed_mean, taed_sd); 

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

file_name = "simulation_dgp1";
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

fid = fopen(bld + "additional_info_dgp1.txt", "w");
fprintf(fid, "%s\n", additional_info{:});
fclose(fid);
