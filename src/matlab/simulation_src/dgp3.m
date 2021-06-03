function [Y, X] = dgp3 (T, N, beta)
    ERROR_SD = sqrt(0.5);
    sigma_sqrd = 1 + 2 *  rand(N * T, 1);
    e          = normrnd(0, ERROR_SD, [T * N, 1]);
    beta       = repmat(beta, [N, 1]);
    
    [X, alpha] = make_X(T, N); 
    Y          = make_Y(X, beta, alpha, sigma_sqrd, e, [], []);
end