function [Y, X] = dgp2 (T, N, beta)
    theta      = repelem(1 + rand(N, 1), T);    
    e          = normrnd(0, sqrt(.75), [T * N, 1]);
    beta       = repmat(beta, [N, 1]);

    [X, alpha] = make_X(T, N); 
    Y          = make_Y(X, beta, alpha, theta, e); 
end 