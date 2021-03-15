function [Y, X1, X2, tau1, tau2, beta1, beta2] = dgp1_autocorr (T, N, a11, a12, a21, a22, sderr_e)
    % returns data and outcomes (X1, X2 and Y) for DGP 1
    
    [beta1, tau1] = make_beta_dgp1(T, 0);
    [beta2, tau2] = make_beta_dgp1(T, 1);
    beta11        = repmat(beta1, [N, 1]);
    beta22        = repmat(beta2, [N, 1]);
    
    alpha      = normrnd(0, 1, [N, 1]);
    alpha      = repelem(alpha, T);
    
    X1 = a11 * normrnd(0, 1, [N * T, 1]) + a12 * alpha;
    X2 = a21 * normrnd(0, 1, [N * T, 1]) + a22 * alpha;

    burn = 50; 
    zeta = normrnd(0, sqrt(sderr_e), [(burn + T), N]); 
    rho  = (ones(N, 1) * .5)'; 
    e    = zeros((T + burn), N);
    
    e(1, 1:N) = zeta(1, 1:N);
    for t = 2:(T + burn)
        e(t, 1:N) = rho .* e(t - 1, 1:N) + zeta(t, 1:N);
    end
	
    e          = e((burn + 1):(burn + T),:);
    e          = reshape(e, T * N, 1);
   
    Y          = alpha + X1 .* beta11 + X2 .* beta22 + e; 
end