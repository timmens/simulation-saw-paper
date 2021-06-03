function [Y, X] = dgp4 (T, N, beta, sd)
    ERROR_SD = sqrt(3);
    if isempty(sd)
        sd = ERROR_SD;
    end
    burn = 100; 
    zeta = normrnd(0, sd, [(burn + T), N]); 
    rho  = (rand(N, 1) * .5 + .25)'; 
    e    = zeros((T + burn), N);
    
    e(1, 1:N) = zeta(1, 1:N);
    for t = 2:(T + burn)
        e(t, 1:N) = rho .* e(t - 1, 1:N) + zeta(t, 1:N);
    end
	
    e          = e((burn + 1):(burn + T),:);
    e          = reshape(e, T * N, 1);
    beta       = repmat(beta, [N, 1]);
    
    [X, alpha] = make_X(T, N); 
    Y          = make_Y(X, beta, alpha, 1, e, [], []);
end