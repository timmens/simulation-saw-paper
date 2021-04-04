function [Y, X, beta] = dgp6 (T, N)
    % returns data and outcomes (X and Y) for DGP 7
    beta   = transpose(repelem(2, T));
    [Y, X] = dgp6(T, N, beta);
end
