function [Y, X] = dgp7 (T, N)
    % returns data and outcomes (X and Y) for DGP 7
    beta   = repelem(2, T);
    [Y, X] = dgp6(T, N, beta);
end