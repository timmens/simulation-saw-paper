% monte carlo study for DGP 1
tic;
        
test = true;
N    = [30, 60]; % , 120, 300];
T    = 2 .^ [5, 6] + 1; %, 7] + 1;
nSim = 4; %500;
nN   = size(N, 2);
nT   = size(T, 2);

option.nGrid     = 40;
option.maxLambda = 100; 
option.minLambda = 0.0001;

s_est_mean = nan(nT * nN, 1);
s_est_sd = nan(nT * nN, 1);
mse_mean1 = nan(nT * nN, 1);
mse_mean2 = nan(nT * nN, 1);
mse_sd1 = nan(nT * nN, 1);
mse_sd2 = nan(nT * nN, 1);
hd_mean1 = nan(nT * nN, 1);
hd_mean2 = nan(nT * nN, 1);
hd_sd1 = nan(nT * nN, 1);
hd_sd2 = nan(nT * nN, 1);

rng_number = 321;
rng(rng_number) % set random number generator seed
parpool(4) % declare (local) parallel cluster with 4 workers 

nIter = nT * nN;

for t = 1:nT
    t_tmp = T(t);
    for n = 1:nN

        disp(strcat('n = ', num2str(N(n)), '; t = ', num2str(T(t))));
        n_tmp = N(n);
        
        results = nan([10, nSim]);
        s_est_mean_tmp = 0;
        s_est_sd_tmp = 0;
        mse_mean_tmp1 = 0;
        mse_mean_tmp2 = 0;
        mse_sd_tmp1  = 0;
        mse_sd_tmp2  = 0;
        hd_mean1_tmp = 0;
        hd_mean2_tmp = 0;
        hd_sd_tmp1   = 0;
        hd_sd_tmp2   = 0;

        parfor r = 1:nSim
            [Y, X1, X2, tau1, tau2, beta1, beta2] = dgp1(t_tmp, n_tmp, 0.5, 0.5, 0.5, 0.5);

            [tau_est,alpha,~,~,~,~,~] = panelpls(Y, [X1, X2], n_tmp, option, 0);  % 1 for mex
            beta_est  = alpha2beta(alpha, tau_est);
            
            tau_est   = tau_est(2:end-1); % they also report first and last time index
            if ~isempty(tau_est)
                tau_est = tau_est - 1;
            end

            s_est_mean_tmp = s_est_mean_tmp + length(tau_est);
            s_est_sd_tmp = s_est_sd_tmp + length(tau_est)^2;

            % mse_mean
            mse_mean_tmp1 = mse_mean_tmp1 + mean((beta1 - beta_est(:, 1)).^2);
            mse_mean_tmp2 = mse_mean_tmp2 + mean((beta2 - beta_est(:, 2)).^2);
            % mse_sd
            mse_sd_tmp1 = mse_sd_tmp1 + mean((beta1 - beta_est(:, 1)).^2)^2;
            mse_sd_tmp2 = mse_sd_tmp2 + mean((beta2 - beta_est(:, 2)).^2)^2;
                                
            % Hausdorff
            hd_mean1_tmp = hd_mean1_tmp + dist_hausdorff(tau1, tau_est);
            hd_mean2_tmp = hd_mean2_tmp + dist_hausdorff(tau2, tau_est);
            % Hausdorff sd
            hd_sd_tmp1 = hd_sd_tmp1 + dist_hausdorff(tau1, tau_est)^2;
            hd_sd_tmp2 = hd_sd_tmp2 + dist_hausdorff(tau2, tau_est)^2;
        end 

        index         = (t - 1) * nN + n;
        s_est_mean(index)  = s_est_mean_tmp / nSim;
        s_est_sd(index) = sqrt((s_est_sd_tmp - s_est_mean_tmp^2 / nSim) / (nSim - 1));
        mse_mean1(index) = mse_mean_tmp1 / nSim;
        mse_mean2(index) = mse_mean_tmp2 / nSim;
        mse_sd1(index) = sqrt((mse_sd_tmp1 - mse_mean_tmp1^2 / nSim) / (nSim - 1));
        mse_sd2(index) = sqrt((mse_sd_tmp2 - mse_mean_tmp2^2 / nSim) / (nSim - 1));
        hd_mean1(index) = hd_mean1_tmp / nSim;
        hd_mean2(index) = hd_mean2_tmp / nSim;
        hd_sd1(index) = sqrt((hd_sd_tmp1 - hd_mean1_tmp^2 / nSim) / (nSim - 1));
        hd_sd2(index) = sqrt((hd_sd_tmp2 - hd_mean2_tmp^2 / nSim) / (nSim - 1));

        % where are we: 
        disp(strcat(num2str((index / nIter) * 100), ' % done'));
    end
end


[cN, cT] = ndgrid(N, T);
T = cT(:); N = cN(:);
result_table = table(T, N, s_est_mean, s_est_sd, mse_mean1, mse_sd1, mse_mean2, mse_sd2, hd_mean1, hd_sd1, hd_mean2, hd_sd2); 

% getting bld path (matlab is such a *** language)
current_file_path = matlab.desktop.editor.getActiveFilename;
bld = split(current_file_path, "/src/");
bld = string(bld);
bld = bld(1);
bld = bld + "/bld/matlab/";

file_name = "simulation_dgp1";
if test
    file_name = file_name + "_" + regexprep(regexprep(datestr(datetime), ' ','-'), ':', '-')+ ".csv";
else
    file_name = file_name + ".csv";
end
writetable(result_table, bld + file_name)

% save additional information
ellapsed_time = toc; 
additional_info        = string(nan(4, 1));
additional_info(1, 1)  = "nsim = " + nSim;
additional_info(2, 1)  = "rng = " + rng_number;
additional_info(3, 1)  = "n.grid = " + option.nGrid;
additional_info(4, 1)  = "ellapsed time = " + ellapsed_time;

fid = fopen(bld + "additional_info_dgp1.txt", "w");
fprintf(fid, "%s\n", additional_info{:});
fclose(fid);
