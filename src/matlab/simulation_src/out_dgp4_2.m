function [Y, X] = dgp4_2 (T, N, beta, var_zeta)     
    % returns data and outcomes (X and Y) for DGP 4 

    % construct autocorrelated error terms (on the individual level);
    % we use a "burn in" which is common in the simulation of stationary
    % time series
    burn = 50; 
    zeta = normrnd(0, sqrt(var_zeta), [(burn + T), N]); 
    rho  = (rand(N, 1) * .5)'; 
    e    = zeros((T + burn), N);
    
    e(1, 1:N) = zeta(1, 1:N);
    for t = 2:(T + burn)
        e(t, 1:N) = rho .* e(t - 1, 1:N) + zeta(t, 1:N);
    end
	
    e          = e((burn + 1):(burn + T),:);
    e          = reshape(e, T * N, 1);
    beta       = repmat(beta, [N, 1]);
    
    [X, alpha] = make_X(T, N); 
    Y          = make_Y(X, beta, alpha, 1, e);
end