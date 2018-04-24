function [vol_loss,X_IBDF,reynolds,u_max,uuu] = Timestepping_SurfTension(L,N,mu,Kernel,X,V,NX,NV,u,tend,dt)

Nx=N(1); Ny=N(2); Nz=N(3);
h = L/Nx;
rho = 1;

vol0 = TetrahedronSphereVolume(X,V,NX,NV);

% set Kernel
kID=Kernel{1}; kID_d=Kernel{2}; K=Kernel{3};

% 3D discrete Laplacian in Fourier space
[m2,m1,m3] = meshgrid(fftshift(-Ny/2:Ny/2-1),...
                      fftshift(-Nx/2:Nx/2-1),...
                      fftshift(-Nz/2:Nz/2-1));
L_hat1 = -4/(h*h)*(sin(pi/Nx*m1).^2+sin(pi/Ny*m2).^2+sin(pi/Nz*m3).^2);
L_hat2 = L_hat1; L_hat2(1,1,1) = 1;

% time step
Nt = floor(tend/dt);
X_IBDF = zeros(NX,3,Nt);
uuu = zeros(Nx,Ny,Nz,3,Nt);
uuu(:,:,:,:,1) = u;
X_IBDF(:,:,1) = X;
vol_loss = zeros(1,Nt);
reynolds = zeros(1,Nt);
u_max    = zeros(1,Nt);

% Nb = 64;
% Ntr = Nb*3;
% th = (0:Nb-1)'*(2*pi)/Nb;
% r1 = 0.125; r2=0.14; r3=0.16;
% Xtracer(1:Nb,:) = [r1*cos(th) + L/2, r1*sin(th)+L/2, zeros(Nb,1)+L/2];
% Xtracer((1:Nb)+Nb,:) = [r2*cos(th)  + L/2,  r2*sin(th)+L/2,  zeros(Nb,1)+L/2];
% Xtracer((1:Nb)+2*Nb,:) = [r3*cos(th)  + L/2,  r3*sin(th)+L/2,  zeros(Nb,1)+L/2];
% load('../sphere/sphere6.mat');
% Xtracer = x(1:nv,:)*0.15 + L/2;
% Vtracer = v;
% volTracer0 = TetrahedronSphereVolume(Xtracer,Vtracer,nv,nt);



        n=0;
        minL = 0.25;
        maxL = 0.75;
        subplot(121)
        th1=trisurf(V,X(:,1),X(:,2),X(:,3), 'LineStyle','none','FaceLighting','phong','FaceAlpha',.6); hold on;
%        th2=trisurf(Vtracer,Xtracer(:,1),Xtracer(:,2),Xtracer(:,3), 'LineStyle','none','FaceLighting','phong','FaceAlpha',.2); hold on;
        set(th1,'Facecolor',[0,0.4470,0.7410]);
%        set(th2,'Facecolor',[0.5,0.5,0.5]);
%         plot3(Xtracer(1:Nb,1),Xtracer(1:Nb,2),Xtracer(1:Nb,3),'.','Markersize',8)
%         plot3(Xtracer((1:Nb)+Nb,1),Xtracer((1:Nb)+Nb,2),Xtracer((1:Nb)+Nb,3),'.','Markersize',8)
%         plot3(Xtracer((1:Nb)+2*Nb,1),Xtracer((1:Nb)+2*Nb,2),Xtracer((1:Nb)+2*Nb,3),'.','Markersize',8)

        light('Position',[0 -2 1])
        colormap('gray')
        axis equal 
        axis([minL,maxL,minL,maxL,minL,maxL])
        xlabel('x')
        ylabel('y')
        zlabel('z')
        view(0,90)
        grid off
        axis on
        hold off
        
        subplot(122)
        th1=trisurf(V,X(:,1),X(:,2),X(:,3), 'LineStyle','none','FaceLighting','phong','FaceAlpha',.6); hold on;
%        th2=trisurf(Vtracer,Xtracer(:,1),Xtracer(:,2),Xtracer(:,3), 'LineStyle','none','FaceLighting','phong','FaceAlpha',.2); hold on;
        set(th1,'Facecolor',[0,0.4470,0.7410]);
%        set(th2,'Facecolor',[0.5,0.5,0.5]);
%         plot3(Xtracer(1:Nb,1),Xtracer(1:Nb,2),Xtracer(1:Nb,3),'.','Markersize',8)
%         plot3(Xtracer((1:Nb)+Nb,1),Xtracer((1:Nb)+Nb,2),Xtracer((1:Nb)+Nb,3),'.','Markersize',8)
%         plot3(Xtracer((1:Nb)+2*Nb,1),Xtracer((1:Nb)+2*Nb,2),Xtracer((1:Nb)+2*Nb,3),'.','Markersize',8)        
        light('Position',[0 -2 1])
        colormap('gray')
        axis equal 
        axis([minL,maxL,minL,maxL,minL,maxL])
        xlabel('x')
        ylabel('y')
        zlabel('z')
        view(55,45)
        grid off
        axis on
        hold off;
        % suptitle(['DFIB, t = ',num2str(n*dt),', Re = ', num2str(reynolds(n))])
        
        
        drawnow;
        
        
        %print('-dpng', sprintf('./figs/surftension_n%d.png', n), '-r300')
        %print('-depsc2', sprintf('./figs/surftension_n%d.eps', n), '-loose')



% use 2nd order Runge-Kutta to get one extra initial condition
% compute X(n+1/2)
u0(1)=mean(mean(mean(u(:,:,:,1)))); 
u0(2)=mean(mean(mean(u(:,:,:,2)))); 
u0(3)=mean(mean(mean(u(:,:,:,3)))); 
Cu = Curl3Dface2edge(u,N,h);
a(:,:,:,1)  = PoissonSolver3D(Cu(:,:,:,1),L_hat2);
a(:,:,:,2)  = PoissonSolver3D(Cu(:,:,:,2),L_hat2);
a(:,:,:,3)  = PoissonSolver3D(Cu(:,:,:,3),L_hat2);
XX = X + (dt/2) * (interpMAC3Dnew(a,X,N,h,kID,kID_d,K)+repmat(u0,NX,1));

% XXtracer = Xtracer + (dt/2) * (interpMAC3Dnew(a,Xtracer,N,h,kID,kID_d,K)+repmat(u0,nv,1));


%force spreading f(n+1/2)
FF = ForceSurfTension(XX,V,NX,NV);
ff = spreadMAC3Dnew(FF,XX,N,h,kID,kID_d,K,L_hat2);
% solve for uu = u(n+1/2)
ADV = advection3D(u,N,h);
c_pressure  = dt/(2*rho);
c_advection = dt/2;
c_Laplacian = 0;
c_force     = dt/(2*rho);
uu = NavierStokes3D_FFT(u,N,h,dt,mu,rho,ADV,ff,L_hat1,L_hat2,c_pressure, c_advection, c_Laplacian, c_force);
% solve for u(n+1)
ADVold = ADV;
ADV = advection3D(uu,N,h);
c_pressure = dt/rho;
c_advection = dt;
c_Laplacian = mu*dt/(2*rho);
c_force = dt/rho;
u = NavierStokes3D_FFT(u,N,h,dt,mu,rho,ADV,ff,L_hat1,L_hat2,c_pressure, c_advection, c_Laplacian, c_force);
ADV = advection3D(u,N,h);
% compute X(n+1)
u0(1)=mean(mean(mean(uu(:,:,:,1)))); 
u0(2)=mean(mean(mean(uu(:,:,:,2)))); 
u0(3)=mean(mean(mean(uu(:,:,:,3))));
Cu = Curl3Dface2edge(uu,N,h);
a(:,:,:,1)  = PoissonSolver3D(Cu(:,:,:,1),L_hat2);
a(:,:,:,2)  = PoissonSolver3D(Cu(:,:,:,2),L_hat2);
a(:,:,:,3)  = PoissonSolver3D(Cu(:,:,:,3),L_hat2);
UU = interpMAC3Dnew(a,XX,N,h,kID,kID_d,K)+repmat(u0,NX,1);
X = X + dt * UU;
X_IBDF(:,:,2) = X;
uuu(:,:,:,:,2) = u;

% Utracer = interpMAC3Dnew(a,Xtracer,N,h,kID,kID_d,K)+repmat(u0,Ntr,1);
% Xtracer = Xtracer + dt * Utracer;
% Xtracer = Xtracer + dt * interpMAC3Dnew(a,XXtracer,N,h,kID,kID_d,K)+repmat(u0,nv,1);

u_max(1)    = max(max(max(sqrt(sum(u.^2,4)))));
vol_loss(1) = abs(TetrahedronSphereVolume(X,V,NX,NV)- vol0) / vol0;
maxUU       = max(sqrt(sum(UU.^2,2)));
reynolds(1) = rho * (L/4) * maxUU / mu;

% volTracer(1) = abs(TetrahedronSphereVolume(Xtracer,Vtracer,nv,nt)- volTracer0) / volTracer0;

        n = 1;
        subplot(121)
        th1=trisurf(V,X(:,1),X(:,2),X(:,3), 'LineStyle','none','FaceLighting','phong','FaceAlpha',.6); hold on;
%         th2=trisurf(Vtracer,Xtracer(:,1),Xtracer(:,2),Xtracer(:,3), 'LineStyle','none','FaceLighting','phong','FaceAlpha',.2); hold on;
        set(th1,'Facecolor',[0,0.4470,0.7410]);
%        set(th2,'Facecolor',[0.5,0.5,0.5]);
%         plot3(Xtracer(1:Nb,1),Xtracer(1:Nb,2),Xtracer(1:Nb,3),'.','Markersize',8)
%         plot3(Xtracer((1:Nb)+Nb,1),Xtracer((1:Nb)+Nb,2),Xtracer((1:Nb)+Nb,3),'.','Markersize',8)
%         plot3(Xtracer((1:Nb)+2*Nb,1),Xtracer((1:Nb)+2*Nb,2),Xtracer((1:Nb)+2*Nb,3),'.','Markersize',8)

        light('Position',[0 -2 1])
        colormap('gray')
        axis equal 
        axis([minL,maxL,minL,maxL,minL,maxL])
        xlabel('x')
        ylabel('y')
        zlabel('z')
        view(0,90)
        grid off
        axis on
        hold off
        
        subplot(122)
        th1=trisurf(V,X(:,1),X(:,2),X(:,3), 'LineStyle','none','FaceLighting','phong','FaceAlpha',.6); hold on;
%        th2=trisurf(Vtracer,Xtracer(:,1),Xtracer(:,2),Xtracer(:,3), 'LineStyle','none','FaceLighting','phong','FaceAlpha',.2); hold on;
        set(th1,'Facecolor',[0,0.4470,0.7410]);
%        set(th2,'Facecolor',[0.5,0.5,0.5]);
%         plot3(Xtracer(1:Nb,1),Xtracer(1:Nb,2),Xtracer(1:Nb,3),'.','Markersize',8)
%         plot3(Xtracer((1:Nb)+Nb,1),Xtracer((1:Nb)+Nb,2),Xtracer((1:Nb)+Nb,3),'.','Markersize',8)
%         plot3(Xtracer((1:Nb)+2*Nb,1),Xtracer((1:Nb)+2*Nb,2),Xtracer((1:Nb)+2*Nb,3),'.','Markersize',8)        
        light('Position',[0 -2 1])
        colormap('gray')
        axis equal 
        axis([minL,maxL,minL,maxL,minL,maxL])
        xlabel('x')
        ylabel('y')
        zlabel('z')
        view(55,45)
        grid off
        axis on
        hold off;
        % suptitle(['DFIB, t = ',num2str(n*dt),', Re = ', num2str(reynolds(n))])
        
        %tightfig
        drawnow;
        
        
        %print('-dpng', sprintf('./figs/surftension_n%d.png', n), '-r300')

for n = 2 : Nt
    
    uold = u;
    
    % compute X(n+1/2)
    u0(1)=mean(mean(mean(u(:,:,:,1)))); 
    u0(2)=mean(mean(mean(u(:,:,:,2)))); 
    u0(3)=mean(mean(mean(u(:,:,:,3))));
    Cu = Curl3Dface2edge(u,N,h);
    a(:,:,:,1)  = PoissonSolver3D(Cu(:,:,:,1),L_hat2);
    a(:,:,:,2)  = PoissonSolver3D(Cu(:,:,:,2),L_hat2);
    a(:,:,:,3)  = PoissonSolver3D(Cu(:,:,:,3),L_hat2);
    XX = X + (dt/2) * (interpMAC3Dnew(a,X,N,h,kID,kID_d,K)+repmat(u0,NX,1));

%    XXtracer = Xtracer + (dt/2) * (interpMAC3Dnew(a,Xtracer,N,h,kID,kID_d,K)+repmat(u0,nv,1));

    
    % compute f(n+1/2)
    FF = ForceSurfTension(XX,V,NX,NV);
    ff = spreadMAC3Dnew(FF,XX,N,h,kID,kID_d,K,L_hat2);

    
    % solve fluid for u(n+1);
    S = 3/2*ADV - 1/2*ADVold;
    c_pressure  = dt/rho;
    c_advection = dt;
    c_Laplacian = mu*dt/(2*rho);
    c_force     = dt/rho;
    u = NavierStokes3D_FFT(u,N,h,dt,mu,rho,S,ff,L_hat1,L_hat2,c_pressure, c_advection, c_Laplacian, c_force);
    
    ADVold = ADV;
    ADV=advection3D(u,N,h);
    
    % compute X(n+1)
    uu = (u+uold)/2;
    uu0(1)=mean(mean(mean(uu(:,:,:,1)))); 
    uu0(2)=mean(mean(mean(uu(:,:,:,2)))); 
    uu0(3)=mean(mean(mean(uu(:,:,:,3))));
    Cuu = Curl3Dface2edge(uu,N,h);
    aa(:,:,:,1)  = PoissonSolver3D(Cuu(:,:,:,1),L_hat2);
    aa(:,:,:,2)  = PoissonSolver3D(Cuu(:,:,:,2),L_hat2);
    aa(:,:,:,3)  = PoissonSolver3D(Cuu(:,:,:,3),L_hat2);
    UU = interpMAC3Dnew(aa,XX,N,h,kID,kID_d,K)+repmat(uu0,NX,1);
    X = X + dt * UU;
    X_IBDF(:,:,n+1) = X;
    uuu(:,:,:,:,n+1) = u;
    
%     Utracer = interpMAC3Dnew(a,Xtracer,N,h,kID,kID_d,K)+repmat(u0,Ntr,1);
%     Xtracer = Xtracer + dt * Utracer;
%    Xtracer = Xtracer + dt * interpMAC3Dnew(aa,XXtracer,N,h,kID,kID_d,K)+repmat(uu0,nv,1);
    
    u_max(n)    = max(max(max(sqrt(sum(u.^2,4)))));
    vol_loss(n) = abs(TetrahedronSphereVolume(X,V,NX,NV)- vol0) / vol0;
    maxUU       = max(sqrt(sum(UU.^2,2)));
    reynolds(n) = rho * (L/4) * maxUU / mu;
    
%    if mod(n,1) == 0
    if ismember(n, [4:4:Nt])
        subplot(121)
        th1=trisurf(V,X(:,1),X(:,2),X(:,3), 'LineStyle','none','FaceLighting','phong','FaceAlpha',.6); hold on;
%        th2=trisurf(Vtracer,Xtracer(:,1),Xtracer(:,2),Xtracer(:,3), 'LineStyle','none','FaceLighting','phong','FaceAlpha',.2); hold on;
        set(th1,'Facecolor',[0,0.4470,0.7410]);
%        set(th2,'Facecolor',[0.5,0.5,0.5]);
%         plot3(Xtracer(1:Nb,1),Xtracer(1:Nb,2),Xtracer(1:Nb,3),'.','Markersize',8)
%         plot3(Xtracer((1:Nb)+Nb,1),Xtracer((1:Nb)+Nb,2),Xtracer((1:Nb)+Nb,3),'.','Markersize',8)
%         plot3(Xtracer((1:Nb)+2*Nb,1),Xtracer((1:Nb)+2*Nb,2),Xtracer((1:Nb)+2*Nb,3),'.','Markersize',8)

        light('Position',[0 -2 1])
        colormap('gray')
        axis equal 
        axis([minL,maxL,minL,maxL,minL,maxL])
        xlabel('x')
        ylabel('y')
        zlabel('z')
        view(0,90)
        grid off
        axis on
        hold off

        subplot(122)
        th1=trisurf(V,X(:,1),X(:,2),X(:,3), 'LineStyle','none','FaceLighting','phong','FaceAlpha',.6); hold on;
%        th2=trisurf(Vtracer,Xtracer(:,1),Xtracer(:,2),Xtracer(:,3), 'LineStyle','none','FaceLighting','phong','FaceAlpha',.2); hold on;
        set(th1,'Facecolor',[0,0.4470,0.7410]);
%        set(th2,'Facecolor',[0.5,0.5,0.5]);
%         plot3(Xtracer(1:Nb,1),Xtracer(1:Nb,2),Xtracer(1:Nb,3),'.','Markersize',10)
%         plot3(Xtracer((1:Nb)+Nb,1),Xtracer((1:Nb)+Nb,2),Xtracer((1:Nb)+Nb,3),'.','Markersize',10)
%         plot3(Xtracer((1:Nb)+2*Nb,1),Xtracer((1:Nb)+2*Nb,2),Xtracer((1:Nb)+2*Nb,3),'.','Markersize',10)        
        light('Position',[0 -2 1])
        axis equal 
        axis([minL,maxL,minL,maxL,minL,maxL])
        xlabel('x')
        ylabel('y')
        zlabel('z')
        view(55,45)
        grid off
        axis on
        hold off;
        %suptitle(['DFIB, t = ',num2str(n*dt),', Re = ', num2str(reynolds(n))])

%        tightfig
        drawnow;

        %print('-depsc2', sprintf('./figs/surftension_n%d.eps', n), '-loose')
        %print('-dpng', sprintf('./figs/surftension_n%d.png', n), '-r300')


    end
    
    
    vol_loss(n) = abs(TetrahedronSphereVolume(X,V,NX,NV) - vol0) / vol0;
%    volTracer(n) = abs(TetrahedronSphereVolume(Xtracer,Vtracer,nv,nt) - volTracer0) / volTracer0;
%    fprintf('n = %d, vol_loss = %.4e, vol_tracer = %.4e, Re = %f, max_u = %f \n', n, vol_loss(n), volTracer(n), reynolds(n), u_max(n))
     fprintf('n = %d, vol_loss = %.4e, Re = %f, max_u = %f \n', n, vol_loss(n), reynolds(n), u_max(n))

end


