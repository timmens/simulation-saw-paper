function [tau] = make_tau_dgp1 (T, plus)
    tau = zeros(2, 1);
    values = [1, 3] + plus;
    for j = 1:2
        tau(j) = floor((values(j) * T) / 5);
    end 
end