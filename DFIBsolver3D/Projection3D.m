% ----------------------------------------------------------------------- %
% Author: Yuanxun Bill Bao
% Date: July 2013
% Description: 3D projection operator
% project a velocity field to its div-free subspace
% i.e., u = Pv = [I-G(DG)^(-1)D]v
%
% Inputs:
% v - Vector fluid velocity 
% N - Domain size [Nx,Ny,Nz].
% L_hat - 3D Laplacian in Fourier space
%
% ----------------------------------------------------------------------- %


function u = Projection3D(v,N,h,L_hat)

Nx = N(1); Ny = N(2); Nz = N(3);

xp = [2:Nx,1]; xm = [Nx,1:Nx-1];
yp = [2:Ny,1]; ym = [Ny,1:Ny-1];
zp = [2:Nz,1]; zm = [Nz,1:Nz-1];

u = zeros(Nx,Ny,Nz,3);

% compute Div v
Divv     = (v(xp,:,:,1) - v(:,:,:,1) + ...
            v(:,yp,:,2) - v(:,:,:,2) + ...
            v(:,:,zp,3) - v(:,:,:,3))/h;
Divv_hat = fftn(Divv);

% (DivGrad)p = Div v 
p_hat    = Divv_hat ./ L_hat;
p        = real(ifftn(p_hat));

% u = v - Grad p
u(:,:,:,1) = v(:,:,:,1) - (p-p(xm,:,:))/h;
u(:,:,:,2) = v(:,:,:,2) - (p-p(:,ym,:))/h;
u(:,:,:,3) = v(:,:,:,3) - (p-p(:,:,zm))/h;

 
