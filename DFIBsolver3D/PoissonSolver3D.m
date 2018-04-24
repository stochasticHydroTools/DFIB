% ----------------------------------------------------------------------- %
% Title: PoissonSolver3D
% Author: Jason Kaye
% Date: June 2013
%
% Description: Solves Poisson's equation with periodic boundary conditions
% in 3D using Fourier basis
%
% Inputs:
% f - Right hand side to Poisson's equation -\Delta u = f, defined on 3D
% grid
% N - Number of grid points for each coordinate, given in form [Nx,Ny,Nz]
% h - Grid point spacing
% L_hat - Discrete Laplacian operator.
%
% Outputs:
% u - Solution defined on same grid as f
%
% Note: Code to generate discrete Laplacian given by
% [m2,m1,m3] = meshgrid(0:N(2)-1, 0:N(1)-1, 0:N(3)-1);
% L_hat = -4/h^2*(sin(pi/N(1)*m1).^2+sin(pi/N(2)*m2).^2+sin(pi/N(3)*m3).^2);
% L_hat(1,1,1) = 1; % Arbitrary choice to avoid divide by zero
% ----------------------------------------------------------------------- %

function u = PoissonSolver3D(f,L_hat)
   
   % Fourier transform of f
   f_hat = fftn(f);
   
   % Solution in Fourier space
   u_hat = f_hat./L_hat;
   
   % Inverse Fourier transform
   u = -real(ifftn(u_hat));
   
end
