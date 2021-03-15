function out = MDCJ (tau, tau_estimate, S)
    out = 0;
    for j = 1:S
        out = out + min(abs(tau_estimate - tau(j)));
    end
end