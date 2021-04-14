function gamma = beta_to_gamma(beta, theta)
    P = size(beta, 2);
    if P > 1
        gamma = [];
        for p = 1:P
            gamma = [gamma, beta(2:end, p), beta(1:end-1, p)];
        end
    else
        gamma = [beta(2:end), beta(1:end-1)];
    end
        
    if isempty(theta)
        gamma = [gamma, zeros(size(gamma, 1), 1)];
    else
        delta_theta = diff(theta);
        gamma = [gamma, delta_theta];
    end
end