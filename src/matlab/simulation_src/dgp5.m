function [Y, X] = dgp5 (T, N, beta)
    gamma = 1 + rand(N * T, 1);
    e = normrnd(0, 1, [T * N, 1]);
    
    beta = repmat(beta, [N, 1]);
    [X, alpha] = make_X(T, N);
    
    theta = make_time_effect(T);
    
    Y = make_Y(X, beta, alpha, gamma, e, theta, []);
end