function [Y, X] = dgp5 (T, N, beta)
    lambda  = normrnd(0, sqrt(.5), [N, 1]);
    f       = normrnd(0, sqrt(.5), [T, 1]);
    epsilon = normrnd(0, sqrt(.5), [T * N, 1]); 
    alpha   = normrnd(0,  1, [N, 1]);
    mu      = normrnd(0,  1, [T * N, 1]);
    
    lambda  = repelem(lambda, T);
    f       = repmat(f, [N, 1]);
    alpha   = repelem(alpha, T);
    
    e       = lambda .* f + epsilon;
    beta    = repmat(beta, [N, 1]);
    
    X       = 0.3 * alpha + 0.3 * lambda + 0.3 * f + mu;
    Y       = make_Y(X, beta, alpha, 1, e);
end