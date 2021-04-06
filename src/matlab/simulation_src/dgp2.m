function [Y, X, Z] = dgp2 (T, N, beta)   

    [~, alpha] = make_X(T, N);
    
    e = normrnd(0, 1, [T * N, 2]);
    variance = [[1, 1/2]; [1/2, 1]];
    L = chol(variance);
    e = e * transpose(L);
    
    Z = alpha / 2 + normrnd(0, 0.1, [N * T, 1]);
    X = 2 * Z + e(:, 1);
    
    beta = repmat(beta, [N, 1]);
    
    Y = make_Y(X, beta, alpha, 1, e(:, 2), [], []);
    
    ZZ = reshape(Z, [T, N]);
    ZZ = ZZ(2:end, :);
    Z = reshape(ZZ, [N * (T - 1), 1]);
end 