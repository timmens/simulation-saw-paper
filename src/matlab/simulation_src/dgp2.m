function [Y, X, Z] = dgp2 (T, N, beta)   

    ERROR_SD = sqrt(0.5);

    [~, alpha] = make_X(T, N);
    
    e = normrnd(0, ERROR_SD, [T * N, 1]);
    
    Z = alpha / 2 + normrnd(0, 1, [N * T, 1]);
    X = 3 * Z + e;
    
    beta = repmat(beta, [N, 1]);
    
    Y = make_Y(X, beta, alpha, 1, e, [], []);
    
    Z0 = reshape(Z, [T, N]);
    Z = [reshape(Z0(2:T,:)-Z0(1:T-1,:),N*(T-1),1) reshape(Z0(2:T,:),N*(T-1),1)];
end
