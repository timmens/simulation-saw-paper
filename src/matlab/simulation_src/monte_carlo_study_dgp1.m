% monte carlo study for DGP 1
profile on;
tic;
        
N    = [25, 50, 100, 200];
T    = 2 .^ [5, 6, 7] + 1;
nSim = 4; %500;
nN   = size(N, 2);
nT   = size(T, 2);

option.nGrid     = 40; 
option.maxLambda = 100; 
option.minLambda = 0.0001;

s_est   = nan(nT * nN, 1);
mse1    = nan(nT * nN, 1);
mse2    = nan(nT * nN, 1);
hd1      = nan(nT * nN, 1);
hd2      = nan(nT * nN, 1);

rng_number = 321;
rng(rng_number)   % set random number generator seed
parpool(4) % declare (local) parallel cluster with 4 workers 

nIter = nT * nN;

for t = 1:nT
    t_tmp = T(t);
    for n = 1:nN
        disp("n = " + N(n) + "; t = " + T(t))
        n_tmp = N(n);
        
        results = nan([10, nSim]);
        s_est_tmp = 0;
        mse_tmp1  = 0;
        mse_tmp2  = 0;
        hd1_tmp    = 0;
        hd2_tmp    = 0;
        for r = 1:nSim
        % parfor r = 1:nSim
            [Y, X1, X2, tau1, tau2, beta1, beta2] = dgp1(t_tmp, n_tmp, 0.5, 0.5, 0.5, 0.5, sqrt(0.75));

            [tau_est,alpha,~,~,~,~,~] = panelpls(Y, [X1, X2], n_tmp, option, 0); % 1 for mex
            beta_est = alpha2beta(alpha, tau_est);
            
            tau_est = tau_est(2:end-1); % they also report first and last time index
            if ~isempty(tau_est)
                tau_est = tau_est - 1;
            end 

            s_est_tmp = s_est_tmp + length(tau_est);

            % MSE
            mse_tmp1   = mse_tmp1  + mean((beta1 - beta_est(:, 1)).^2);
            mse_tmp2   = mse_tmp2  + mean((beta2 - beta_est(:, 2)).^2);
            
                                
            % Hausdorff
            hd1_tmp    = hd1_tmp + dist_hausdorff(tau1, tau_est);
            hd2_tmp    = hd2_tmp + dist_hausdorff(tau2, tau_est);
        end 

        index         = (t - 1) * nN + n;
        s_est(index)  = s_est_tmp / nSim;
        mse1(index)   = mse_tmp1 / nSim;
        mse2(index)   = mse_tmp2 / nSim;
        hd1(index)    = hd1_tmp / nSim;
        hd2(index)    = hd2_tmp / nSim;

        % where are we: 
        disp((index / nIter) * 100 + " % done");
    end
end

% format table and write to disc
[cN, cT] = ndgrid(N, T);
T = cT(:); N = cN(:);
result_table = table(T, N, s_est, mse1, mse2, hd1, hd2); 
disp(result_table);
ellapsed_time = toc;

% getting bld path (matlab is such a *** language)
current_file_path = matlab.desktop.editor.getActiveFilename;
bld = split(current_file_path, "/src/");
bld = string(bld);
bld = bld(1);
bld = bld + "/bld/matlab/";

file_name = "simulation-dgp1-" + regexprep(regexprep(datestr(datetime), ' ','-'), ':', '-')+ ".csv";
writetable(result_table, bld + file_name)

% save additional information
additional_info        = string(nan(4, 1));
additional_info(1, 1)  = "nsim = " + nSim;
additional_info(2, 1)  = "rng = " + rng_number;
additional_info(3, 1)  = "n.grid = " + option.nGrid;
additional_info(4, 1)  = "ellapsed time = " + ellapsed_time;

fid = fopen(bld + "additional_info_dgp1.txt", "w");
fprintf(fid, "%s\n", additional_info{:});
fclose(fid);

profile viewer