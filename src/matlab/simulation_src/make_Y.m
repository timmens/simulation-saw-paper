function out = make_Y (X, beta, alpha, gamma, e, theta, mu)
    if isempty(mu)
        mu = 0;
    end
    if isempty(theta)
        theta = 0;
    else
        N = length(alpha) / length(theta);
        theta = repmat(theta, N, 1);
    end
    
	out = mu + X .* beta + alpha + theta + sqrt(gamma) .* e;
end