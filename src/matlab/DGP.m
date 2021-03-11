function [Y, X] = DGP (T, N, beta, index)
    if index <= 0
        error("index has to be a whole number between 2 and 6");
    elseif index == 1
        error("dgp 1 cannot be called through this wrapper function");
    elseif index == 2
        [Y, X] = dgp2(T, N, beta);
    elseif index == 3
        [Y, X] = dgp3(T, N, beta);
    elseif index == 4
        [Y, X] = dgp4(T, N, beta);
    elseif index == 5
        [Y, X] = dgp5(T, N, beta);    
    elseif index == 6
        [Y, X] = dgp6(T, N, beta);    
    elseif index == 7
        [Y, X] = dgp7(T, N); %error("dgp 7 cannot be called through this wrapper function");
    else
        error("index has to be a whole number between 2 and 6");
    end