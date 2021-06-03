function [beta, tau] = make_beta (T, S, N)
    % constructs the beta jump function using the number of 
    % time periods and the number of total jumps
    
    % returns beta vector and tau (jump position set)
    
    if N == 30
        magnitude = 7;
    elseif N == 60
        magnitude = 5;
    elseif N == 120
        magnitude = 4;
    elseif N == 300
        magnitude = 3;
    end
    
    magnitude = magnitude / 3;
    
    if S == 0
       beta = repelem(magnitude, T)'; 
    else 
       tau      = make_tau(T, S);
       rep_beta = diff([0; tau; T]);
       betas    = magnitude * (-1).^(1:(S + 1)).';  
       beta     = repelem(betas, rep_beta);
    end
    
end
