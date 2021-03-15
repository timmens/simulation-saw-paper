function [beta, tau] = make_beta_dgp1 (T, plus)
    % constructs the beta jump function using the number of 
    % time periods and the number of total jumps
    
    % returns beta vector and tau (jump position set)
    tau      = make_tau_dgp1(T, plus);
    rep_beta = diff([0; tau; T]);
    betas    = (2 / 3) * (-1).^(1:(2 + 1)).';  
    
    beta      = repelem(betas, rep_beta);
    
end