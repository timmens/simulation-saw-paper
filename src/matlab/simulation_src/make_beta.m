function [beta, tau] = make_beta (T, S)
    % constructs the beta jump function using the number of 
    % time periods and the number of total jumps
    
    % returns beta vector and tau (jump position set)
    
    if S == 0
       beta = repelem(2, T)'; 
    else 
       tau      = make_tau(T, S);
       rep_beta = diff([0; tau; T]);
       betas    = 2 * (-1).^(1:(S + 1)).';  
       beta      = repelem(betas, rep_beta);
    end
    
end
