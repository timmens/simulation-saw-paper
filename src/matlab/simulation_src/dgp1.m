function [Y, X1, X2, tau1, tau2, beta1, beta2] = dgp1 (T, N, a11, a12, a21, a22, sderr_e)
    % returns data and outcomes (X1, X2 and Y) for DGP 1
    
    [beta1, tau1] = make_beta_dgp1(T, 0);
    [beta2, tau2] = make_beta_dgp1(T, 1);
    beta11        = repmat(beta1, [N, 1]);
    beta22        = repmat(beta2, [N, 1]);
    
    alpha      = normrnd(0, 1, [N, 1]);
    alpha      = repelem(alpha, T);
    
    X1 = a11 * normrnd(0, 1, [N * T, 1]) + a12 * alpha;
    X2 = a21 * normrnd(0, 1, [N * T, 1]) + a22 * alpha;
    
    e          = normrnd(0, sderr_e, [T * N, 1]);
   
    Y          = alpha + X1 .* beta11 + X2 .* beta22 + e; 
end