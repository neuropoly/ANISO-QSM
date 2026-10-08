function [ output ] = gradient( m )
    N = size(m);
    [k1, k2, k3] = ndgrid(0:N(1)-1, 0:N(2)-1, 0:N(3)-1);
    E1 = 1 - exp(2i .* pi .* k1 / N(1));
    E2 = 1 - exp(2i .* pi .* k2 / N(2));
    E3 = 1 - exp(2i .* pi .* k3 / N(3)); 
    
    % Convolution
    Gm1 = real(ifftn(E1.*fftn(m)));
    Gm2 = real(ifftn(E2.*fftn(m)));
    Gm3 = real(ifftn(E3.*fftn(m)));
    
    % Generate output
    output = abs( Gm1 ) + abs( Gm2 ) + abs( Gm3 );
    
end

