function out = make_Y (X, beta, alpha, theta, e)
	out = X .* beta + alpha + sqrt(theta) .* e;
end