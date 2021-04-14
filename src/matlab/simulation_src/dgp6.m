function [Y, X, beta] = dgp6 (T, N)
    beta   = transpose(repelem(2, T));
    [Y, X] = dgp4(T, N, beta);
end
