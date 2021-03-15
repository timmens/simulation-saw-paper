function [X, alpha] = make_X (T, N) 
    % constructs the data matrix;
    % individual effect (alpha) is also appearing in 
    % the outcome, hence they are returned as well

    xi    = normrnd(0, 1, [T * N, 1]);
    alpha = normrnd(0, 1, [N, 1]);
    alpha = repelem(alpha, T);
    
    X     = 0.5 * alpha + xi;
end