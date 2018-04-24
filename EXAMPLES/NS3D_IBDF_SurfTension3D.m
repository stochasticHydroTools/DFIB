clear; close all; restoredefaultpath;
addpath('../Kernels_MEX',...
        '../DFIB_SpreadInterp3D_MEX',...
        '../DFIBsolver3D');

sn = 'sphere4';
load( sprintf('%s.mat',sn) )

profile off;                                    

L = 1; 
mu = 0.05;

%Kernel={'bspline4pt', 'bspline4pt_d',[]};

K= (59/60)*(1-sqrt(1-(3220/3481)));
Kernel = {'flex6pt','flex6pt_d',K};

Nx = 64; Ny = Nx; Nz=Nx;
N = [Nx,Ny,Nz];
h = L/Nx;

MpC = 1;
hIB = h/MpC;

R = sqrt(sqrt(3)/4 *hIB^2 * nt / (4*pi));

% Lagrangian mesh
NX = nv; NV = nt;
X = x(1:NX,:)*R + L/2;
V = v;

% initial velocity
x = (0:Nx-1)*h;
y = (0:Ny-1)*h;
z = (0:Nz-1)*h;
[yy1,xx1,zz1] = meshgrid(y+h/2,x,z+h/2);
[yy2,xx2,zz2] = meshgrid(y,x+h/2,z+h/2);
[yy3,xx3,zz3] = meshgrid(y+h/2,x+h/2,z);
[yy,xx,zz] = meshgrid(y+h/2,x+h/2,z+h/2);

initial = 'shear'; 

v = zeros(Nx,Ny,Nz,3);
switch (initial)
    case 'zero'
        v = zeros(Nx,Ny,Nz,3);
    case 'shear'
        v(:,:,:,1) = 0;
        v(:,:,:,2) = 2*sin(2*pi*xx2/L);
        v(:,:,:,3) = 0;
    case 'taylor'
        p = 2;
        v(:,:,:,1) =  8*cos(2*pi*xx1/L*p)  .* sin(2*pi*yy1/L*p) .* sin(2*pi*zz1/L*p);
        v(:,:,:,2) = -4*sin(2*pi*xx2/L*p) .* cos(2*pi*yy2/L*p) .* sin(2*pi*zz2/L*p);
        v(:,:,:,3) = -4*sin(2*pi*xx3/L*p) .* sin(2*pi*yy3/L*p) .* cos(2*pi*zz3/L*p);
end

% 3D Laplacian in Fourier Space
[m2,m1,m3] = meshgrid(fftshift(-Ny/2:Ny/2-1),...
                      fftshift(-Nx/2:Nx/2-1),...
                      fftshift(-Nz/2:Nz/2-1));
L_hat = -4/(h*h)*(sin(pi/Nx*m1).^2+sin(pi/Ny*m2).^2+sin(pi/Nz*m3).^2);
L_hat(1,1,1) = 1;
% project initial to div-free
ustart = Projection3D(v,N,h,L_hat);

% time step
tend = 0.5;
dt   = h*0.25;
Nt   = floor(tend/dt);
dt   = tend/Nt;
tt   = 0:dt:tend;

%fprintf('# of markers per cell = %f \n', MpC);

tic
[vol_lossIBDF,X_IBDF,Re,u_max,uvw] = Timestepping_SurfTension(L,N,mu,Kernel,X,V,NX,NV,ustart,tend,dt);
% FILEname = sprintf('/scratch/billbao/NS3D-IBDF-%s-%s-MpC%.2f-N%d-%s.mat',sn,initial,MpC,Nx,Kernel{1});
% FILEname = sprintf('~/Desktop/DivIBdata/paraview-%s-%s-MpC%.2f-N%d-%s.mat',sn,initial,MpC,Nx,Kernel{1});
% save(FILEname,'tt','vol_lossIBDF','X_IBDF','uvw','V','h','MpC','Nx','L','mu','R','Re','u_max','xx','yy','zz');

toc

% profile viewer
% p = profile('info');
% profsave(p,'profile_results')
% profile off;
