function [Y, X1, X2, tau1, tau2, beta1, beta2] = dgp1 (T, N)
    % returns data and outcomes (X1, X2 and Y) for DGP 1

    ERROR_SD = sqrt(0.5);
    
    [beta1, tau1] = make_beta(T, 2, N);
    [beta2, tau2] = make_beta(T, 3, N);
    beta11        = repmat(beta1, [N, 1]);
    beta22        = repmat(beta2, [N, 1]);
    
    alpha      = normrnd(0, 1, [N, 1]);
    alpha      = repelem(alpha, T);
    
    X1 = 0.5 * alpha + normrnd(0, 1, [N * T, 1]);
    X2 = 0.5 * alpha + normrnd(0, 1, [N * T, 1]);
    
    e = normrnd(0, ERROR_SD, [T * N, 1]);
   
    Y = alpha + X1 .* beta11 + X2 .* beta22 + e; 
end
