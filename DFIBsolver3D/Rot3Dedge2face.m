%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rot3Dedge2face.m
% July 13, 2013
% Yuanxun Bill Bao
%
% 3D discrete Rot operator: 
% edge-centered vector field -> face-centered vector field
%
% Rot v = [ D2v3 - D3v2,
%           D3v1 - D1v3,
%           D1v2 - D2v1 ]
%
% For stencils, refer to RotCurl.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function u = Rot3Dedge2face(v,N,h)

Nx = N(1); Ny = N(2); Nz = N(3);

xp = [2:Nx,1]; 
yp = [2:Ny,1]; 
zp = [2:Nz,1];

u = zeros(Nx,Ny,Nz,3);

% D2v3 - D3v2
u(:,:,:,1) = (v(:,yp,:,3) - v(:,:,:,3))/h - ...
             (v(:,:,zp,2) - v(:,:,:,2))/h;

% D3v1 - D1v3
u(:,:,:,2) = (v(:,:,zp,1) - v(:,:,:,1))/h - ...
             (v(xp,:,:,3) - v(:,:,:,3))/h; 

% D1v2 - D2v1         
u(:,:,:,3) = (v(xp,:,:,2) - v(:,:,:,2))/h - ...
             (v(:,yp,:,1) - v(:,:,:,1))/h;  

