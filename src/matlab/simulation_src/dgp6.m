function [Y, X] = dgp6 (T, N, beta)
    % returns data and outcomes (X and Y) for DGP 6
    
    burn    = 50; 
    zeta_mu = normrnd(0, sqrt(.5), [(burn + T), N]);
    zeta_e  = normrnd(0, sqrt(.5), [(burn + T), N]); 
    rho_e   = (rand(N, 1) * .5)'; 
    rho_mu  = (rand(N, 1) * .5)';
    mu      = zeros((T + burn), N);
    epsilon = zeros((T + burn), N);
    
    mu(1, 1:N)      = zeta_mu(1, 1:N);
    epsilon(1, 1:N) = zeta_e(1, 1:N);
    for t = 2:(T + burn)
        mu(t, 1:N)      = rho_mu .* mu(t - 1, 1:N)      + zeta_mu(t, 1:N);
        epsilon(t, 1:N) = rho_e  .* epsilon(t - 1, 1:N) + zeta_e(t, 1:N);
    end
	
    mu      = mu((burn + 1):(burn + T),:);
    mu      = reshape(mu, T * N, 1);
    epsilon = epsilon((burn + 1):(burn + T),:);
    epsilon = reshape(epsilon, T * N, 1);
    	
    
    lambda  = normrnd(0, sqrt(.5), [N, 1]);
    f       = normrnd(0, sqrt(.5), [T, 1]);
    alpha   = normrnd(0,  1, [N, 1]);
    
    lambda  = repelem(lambda, T);
    f       = repmat(f, [N, 1]);
    alpha   = repelem(alpha, T);
    
    e       = lambda .* f + epsilon;
    beta    = repmat(beta, [N, 1]);
    
    X       = 0.3 * alpha + 0.3 * lambda + 0.3 * f + mu;
    Y       = make_Y(X, beta, alpha, 1, e);
end