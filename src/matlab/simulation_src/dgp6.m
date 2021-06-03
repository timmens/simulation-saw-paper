function [Y, X, beta] = dgp6 (T, N)
    beta = transpose(repelem(1, T));
    [Y, X] = dgp4(T, N, beta, 2);
end
