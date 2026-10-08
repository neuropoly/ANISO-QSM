function [ kernel ] = dipole_kernel( N, spatial_res)

[ky,kx,kz] = (meshgrid(-floor(N(2)/2):ceil(N(2)/2)-1, -floor(N(1)/2):ceil(N(1)/2)-1, -floor(N(3)/2):ceil(N(3)/2)-1));

kx = (single(kx) / max(abs(single(kx(:))))) / spatial_res(1);
ky = (single(ky) / max(abs(single(ky(:))))) / spatial_res(2);
kz = (single(kz) / max(abs(single(kz(:))))) / spatial_res(3);

k2 = kx.^2 + ky.^2 + kz.^2;
k2(k2==0) = eps;

kernel = ifftshift( 1/3 - (kz.^2) ./ k2 ); 
kernel(1,1,1) = 0.0;

