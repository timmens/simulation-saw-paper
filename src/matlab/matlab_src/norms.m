function cvx_optval = norms( x, p, dim )
switch p,
    case 1,
        cvx_optval = sum( abs( x ), dim );
    case 2,
        cvx_optval = sqrt( sum( x .* conj( x ), dim ) );
    case Inf,
        cvx_optval = max( abs( x ), [], dim );
    otherwise,
        cvx_optval = sum( abs( x ) .^ p, dim ) .^ ( 1 / p );
end
