function [Y, X] = dgp3 (T, N, beta)
    theta      = 1 + rand(N * T, 1);
    e          = normrnd(0, sqrt(.75), [T * N, 1]); % normrnd(0, sqrt(.75), [T * N, 1]); 
    beta       = repmat(beta, [N, 1]);
    
    [X, alpha] = make_X(T, N); 
    Y          = make_Y(X, beta, alpha, theta, e);
end   