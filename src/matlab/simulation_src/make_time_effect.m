function theta = make_time_effect (T)
    [theta, ~] = make_beta(T, floor(T / 10), 30);
end 