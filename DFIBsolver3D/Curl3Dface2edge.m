%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Curl3D2face2edge.m
% July 13, 2013
% Yuanxun Bill Bao
%
% 3D discrete Curl operator: 
% face-centered vector field -> edge-centered vector field
%
% Curl v = [ D2v3 - D3v2,
%            D3v1 - D1v3,
%            D1v2 - D2v1 ]
%
% For stencils, refer to RotCurl.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = Curl3Dface2edge(v,N,h)

Nx = N(1); Ny = N(2); Nz = N(3);

xm = [Nx,1:Nx-1]; 
ym = [Ny,1:Ny-1]; 
zm = [Nz,1:Nz-1];

u = zeros(Nx,Ny,Nz,3);

u(:,:,:,1) = (v(:,:,:,3) - v(:,ym,:,3))/h - ...
             (v(:,:,:,2) - v(:,:,zm,2))/h;

u(:,:,:,2) = (v(:,:,:,1) - v(:,:,zm,1))/h - ...
             (v(:,:,:,3) - v(xm,:,:,3))/h; 
         
u(:,:,:,3) = (v(:,:,:,2) - v(xm,:,:,2))/h - ...
             (v(:,:,:,1) - v(:,ym,:,1))/h;  

