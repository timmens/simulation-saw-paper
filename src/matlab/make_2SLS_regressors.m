% function X_proj = make_2SLS_regressors (X, T, k)
%    N     = length(X) / T;
%    X     = reshape(X, [N, T]);
%    Z     = make_instruments(X, k);
%    
%    X_proj       = zeros(N, T);
%    for t = 1:T
%        X_proj(:, t) = predict(fitlm(Z, X(:, t)));
%    end
%    X_proj = reshape(X_proj, [N * T, 1]);
% end
