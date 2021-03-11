function res = make_aux_iv (X,T,N,k)
    X = reshape(X, [T, N]);
    [gamma, ~] = eig(cov(X));
    if k > size(gamma, 2)
        k = size(gamma, 2);
    end
    gamma_k = gamma(:, 1:k);
    
    X_centered = X - mean(X);
    beta_mat   = X_centered * gamma_k;
    X_cent_hat = beta_mat * transpose(gamma_k);
    
    res = X_centered - X_cent_hat;
    res = reshape(res, [T * N, 1]);
end