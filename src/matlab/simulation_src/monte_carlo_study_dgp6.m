profile on;
tic;
     
N    = [30, 60, 120, 300];
T    = 2 .^ [5, 6, 7] + 1;
nSim = 500;
nN   = size(N, 2);
nT   = size(T, 2);

option.nGrid     = 50; 
option.maxLambda = 100; 
option.minLambda = 0.0001;

s_est = nan(nT * nN, 1);
mdcj  = nan(nT * nN, 1);
mise  = nan(nT * nN, 1);
hd    = nan(nT * nN, 1);
s_0   = nan(nT * nN, 1);

rng_number = 123;
rng(rng_number) % set random number generator seed
parpool(4) % declare (local) parallel cluster with 4 workers 

nIter = nT * nN;

for t = 1:nT
    t_tmp = T(t);
    beta = repmat(2, [T, 1]);
    for n = 1:nN
        n_tmp = N(n);
        disp("n = " + n_tmp + "; s = " + s_tmp + "; t = " + t_tmp)
            
        s_est_tmp = 0;
        mdcj_tmp  = 0;
        mise_tmp  = 0;
        hd_tmp    = 0;
        s_0_tmp   = 0;
        
        parfor r = 1:nSim
            [Y, X] = dgp6(t_tmp, n_tmp);
            
            [tau_est,alpha_est,~,~,~,~,~] = panelpls(Y, X, n_tmp, option, 1);
            
            beta_est  = alpha2beta(alpha_est, tau_est);
            
            tau_est   = tau_est(2:end-1) % they also report first and last time index
            if ~isempty(tau_est)
                tau_est = tau_est - 1
                % they report tau_i as starting of next regime we as end of current regime
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%% CODE EXCEPTIONS %%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % s_est
            s_est_tmp = s_est_tmp + length(tau_est);
            
            % MDCJ
            mdcj_tmp  = mdcj_tmp + MDCJ(tau, tau_est, S);
            
            % MISE
            mise_tmp  = mise_tmp + mean((beta - beta_est).^2);
            
            % Hausdorff
            hd_tmp    = hd_tmp + dist_hausdorff(tau, tau_est);
            
            % Check if no jump location is found
            if isempty(tau_est)
                s_0_tmp = s_0_tmp + 1;
            end
        end
        
        index = (t - 1) * nS * nN + (s - 1) * nN + n;
        s_est(index) = s_est_tmp  / nSim;
        mdcj(index)  = mdcj_tmp  / nSim;
        mise(index)  = mise_tmp / nSim;
        s_0(index)   = s_0_tmp / nSim;
        hd(index)    = hd_tmp / nSim;
        
        % where are we:
        disp((index / nIter) * 100 + " % done");
    end
end

% getting bld path (matlab is such a *** language)
current_file_path = matlab.desktop.editor.getActiveFilename;
bld = split(current_file_path, "/src/");
bld = string(bld);
bld = bld(1);
bld = bld + "/bld/matlab/";

file_name = "simulation-dgp6-" + regexprep(regexprep(datestr(datetime), ' ','-'), ':', '-')+ ".csv";
writetable(result_table, bld + file_name)

% save additional information
additional_info        = string(nan(4, 1));
additional_info(1, 1)  = "nsim = " + nSim;
additional_info(2, 1)  = "rng = " + rng_number;
additional_info(3, 1)  = "n.grid = " + option.nGrid;
additional_info(4, 1)  = "ellapsed time = " + ellapsed_time;

fid = fopen(bld + "additional_info_dgp6.txt", "w");
fprintf(fid, "%s\n", additional_info{:});
fclose(fid);