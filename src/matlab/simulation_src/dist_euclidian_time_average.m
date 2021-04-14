function res = dist_euclidian_time_average (gamma_hat, gamma)
    dist = 0;
    T = size(gamma, 1);
    for t = 1:T
        dist = dist + sum((gamma_hat(t, :) - gamma(t, :)).^2);
    end
    
    res = dist / T;
end