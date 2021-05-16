function out = make_tau (T, S)
    % computes jump positions given the number of time periods 
    % and the number of total jumps, via the algorithm given in
    % the paper 

    out = zeros(S, 1);
    for j = 1:S
        out(j) = floor((T - 1) / (2^j)) + 1;
    end
    out = sort(out);
end
 