% ----------------------------------------------------------------------- %
% NavierStokes3D_FFT
% Yuanxun Bill Bao
% July 1, 2014
%
% Description: FFT based 2D Navier Stokes solver
% 
% Inputs:
% N      -- grid size N=[Nx,Ny,Nz], Nx=Ny=Nz for now
% mu     -- viscosity
% rho    -- density
% S      -- advection term
% f      -- forcing
% L_hat1 -- Laplacian in Fourier space
% L_hat2 -- Laplacian in Fourier space with zeroth mode fixed to be 1
% ----------------------------------------------------------------------- %


function unew = NavierStokes3D_FFT(u,N,h,dt,mu,rho,S,f,L_hat1,L_hat2,...
                c_pressure, c_advection, c_Laplacian, c_force)

Nx = N(1); Ny = N(2); Nz = N(3);

c = mu*dt/(2*rho);

% poissone solve for the pressure p(n+1/2):
w1 = u(:,:,:,1)-c_advection*S(:,:,:,1)+c_Laplacian*Laplacian3D(u(:,:,:,1),N,h)+c_force*f(:,:,:,1);
w2 = u(:,:,:,2)-c_advection*S(:,:,:,2)+c_Laplacian*Laplacian3D(u(:,:,:,2),N,h)+c_force*f(:,:,:,2);
w3 = u(:,:,:,3)-c_advection*S(:,:,:,3)+c_Laplacian*Laplacian3D(u(:,:,:,3),N,h)+c_force*f(:,:,:,3);
Divw = (w1([2:Nx,1],:,:)-w1 + w2(:,[2:Ny,1],:)-w2 + w3(:,:,[2:Nz,1])-w3)/h;
pp = -1/c_pressure * PoissonSolver3D(Divw,L_hat2);

% solve for velocity u(n+1)
r1 = u(:,:,:,1) + c_Laplacian*Laplacian3D(u(:,:,:,1),N,h) - c_advection*S(:,:,:,1) - ...
     c_pressure*(pp-pp([Nx,1:Nx-1],:,:))/h + c_force*f(:,:,:,1);
r2 = u(:,:,:,2) + c_Laplacian*Laplacian3D(u(:,:,:,2),N,h) - c_advection*S(:,:,:,2) - ...
     c_pressure*(pp-pp(:,[Ny,1:Ny-1],:))/h + c_force*f(:,:,:,2);
r3 = u(:,:,:,3) + c_Laplacian*Laplacian3D(u(:,:,:,3),N,h) - c_advection*S(:,:,:,3) - ...
     c_pressure*(pp-pp(:,:,[Nz,1:Nz-1]))/h + c_force*f(:,:,:,3);
u1_hat = fftn(r1) ./ (1-c*L_hat1);
u2_hat = fftn(r2) ./ (1-c*L_hat1);
u3_hat = fftn(r3) ./ (1-c*L_hat1);

unew(:,:,:,1) = real(ifftn(u1_hat));
unew(:,:,:,2) = real(ifftn(u2_hat));
unew(:,:,:,3) = real(ifftn(u3_hat));

