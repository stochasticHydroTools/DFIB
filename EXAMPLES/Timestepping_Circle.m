% ----------------------------------------------------------------------- %
% NS2D_IBDF.m
% Yuanxun Bill Bao
% Nov 09, 2013
% 
% Description: solve the fluid-structure interaction of an elastic
% membrane by solving a 2D NS flow. The new immersed boundary
% method with divergence-free interpolated velocity (IB-DF) is used. 
% A midpoint temporal integrator is used for IB markers
% A second-order AB scheme is used for fluid solve
%
% Inputs:
% L      -- domain length
% N      -- grid size N=[Nx,Ny], Nx=Ny for now
% mu     -- viscosity
% Kernel -- [KernelID, DerivativeID, K]
% X0     -- initial configuration of the memberane
% u0     -- intial velocity
% tend   -- final time
% dt     -- time step
% ----------------------------------------------------------------------- %

function [uIBDF,XIBDF,areaIBDF] = Timestepping_Circle(L,N,mu,rho,Kernel,X,u,tend,dt,Nf,MpC,showplot)
%function [areaIBDF,areaTracer,areaTracer_para,TH,Fpol,Xtracer] = NS2D_IBDF_midpoint(L,N,mu,Kernel,X,u,tend,dt,Nf,MpC,showplot)
%function [areaIBDF,areaTracer,u,X,Xtracer,TH,Fpol] = NS2D_IBDF_midpoint(L,N,mu,Kernel,X,u,tend,dt,Nf,MpC,showplot)

Nx=N(1); Ny=N(2);
h=L/Nx;
Ns=length(X); 
ds=2*pi/Ns;

x=(0:Nx-1)*h;
y=(0:Ny-1)*h;
[yy,xx]=meshgrid(y,x);
xx2=xx+h/2; yy2=yy+h/2;
sk=1:1:Nx;

kp = [2:Nx,    1];
km = [Nx  ,1:Nx-1];
xp = [2:Nx,1];
xm = [Nx,1:Nx-1];
yp = [2:Ny,1];
ym = [Ny,1:Ny-1];

% set Kernel
kID=Kernel{1}; kID_d=Kernel{2}; K=Kernel{3};

% 2D Laplacian in Fourier space
[m2,m1] = meshgrid(fftshift(-Ny/2:Ny/2-1), fftshift(-Nx/2:Nx/2-1));
L_hat1  = -4/(h*h) * (sin(pi/N(1)*m1).^2 + sin(pi/N(2)*m2).^2);
L_hat2  = L_hat1; L_hat2(1,1)=1;

% % add passive tracers
% Ntracer = Ns*20;
% s = (0:Ntracer-1)*(2*pi/Ntracer);
% alpha=1/4; beta=1/4;
% %alpha = 5/28; beta = 7/20;
% Xtracer = [alpha*cos(s'), beta*sin(s')]*L+L/2;
% Xtracer0 = Xtracer;

% sk=1;
% areaTracer(1,1) = polyarea([Xtracer(1:sk:end,1); Xtracer(1,1)], [Xtracer(1:sk:end,2);Xtracer(1,2)]);
% [ampl,ppint,xf,yf] = makeIBcurve(Xtracer(1:sk:end,1)-L/2,Xtracer(1:sk:end,2)-L/2,100);
% areaTracer(2,1) = ppint;
% sk=2;
% areaTracer(3,1) = polyarea([Xtracer(1:sk:end,1); Xtracer(1,1)], [Xtracer(1:sk:end,2);Xtracer(1,2)]);
% [ampl,ppint,xf,yf] = makeIBcurve(Xtracer(1:sk:end,1)-L/2,Xtracer(1:sk:end,2)-L/2,100);
% areaTracer(4,1) = ppint;
% sk=20;
% areaTracer(5,1) = polyarea([Xtracer(1:sk:end,1); Xtracer(1,1)], [Xtracer(1:sk:end,2);Xtracer(1,2)]);
% [ampl,ppint,xf,yf] = makeIBcurve(Xtracer(1:sk:end,1)-L/2,Xtracer(1:sk:end,2)-L/2,100);
% areaTracer(6,1) = ppint;

% time step
Nt = floor(tend/dt);

%uIBDF=zeros(Nx,Ny,2,Nf+1); XIBDF=zeros(Ns,2,Nf+1);
areaIBDF=zeros(1,Nt+1);
uIBDF(:,:,:,1)=u;
XIBDF(:,:,1)=X;
areaIBDF(1,1) = polyarea([X(:,1); X(1,1)], [X(:,2);X(1,2)]);
[ampl,ppint,xf,yf] = makeIBcurve(X(:,1)-L/2,X(:,2)-L/2,Ns);
areaIBDF(2,1) = ppint;

%% use 2nd order Runge-Kutta to get one extra initial value
% compute X(n+1/2)
u0(1)=mean(mean(u(:,:,1))); u0(2)=mean(mean(u(:,:,2)));
Cu = Curl2Dedge2node(u,N,h);
a  = PoissonSolver2D(Cu,L_hat2);
XX = X + (dt/2) * (interpMAC2Dnew(a,X,N,h,kID,kID_d,K)+repmat(u0,Ns,1));
%XXtracer = Xtracer + (dt/2) * (interpMAC2Dnew(a,Xtracer,N,h,kID,kID_d,K)+repmat(u0,Ntracer,1));
% Force spreading f(n+1/2)
FF=Force_HookeSpring(XX,Ns,ds)*ds;
ff = spreadMAC2Dnew(FF,XX,N,h,kID,kID_d,K,L_hat2);
% solve for uu=u(n+1/2)
ADV=advection2D(u,N,h);
w1=(-2*rho/dt)*u(:,:,1)+rho*ADV(:,:,1);%-ff(:,:,1);
w2=(-2*rho/dt)*u(:,:,2)+rho*ADV(:,:,2);%-ff(:,:,2);
Divw=(w1(xp,:)-w1)/h + (w2(:,yp)-w2)/h;
pp=PoissonSolver2D(Divw,L_hat2);
c1=dt/(2*rho); c2=mu*dt/(2*rho);
r1=u(:,:,1)-(dt/2)*ADV(:,:,1)+c1*ff(:,:,1)-c1*(pp-pp(xm,:))/h;
r2=u(:,:,2)-(dt/2)*ADV(:,:,2)+c1*ff(:,:,2)-c1*(pp-pp(:,ym))/h;
uu1_hat=fft2(r1)./(1-c2*L_hat1);
uu2_hat=fft2(r2)./(1-c2*L_hat1);
uu(:,:,1)=real(ifft2(uu1_hat));
uu(:,:,2)=real(ifft2(uu2_hat));

% solve for u(n+1)
ADVold = ADV;
ADV=advection2D(uu,N,h);
w1=(-rho/dt)*u(:,:,1)+rho*ADV(:,:,1)-mu/2*Laplacian2D(u(:,:,1),N,h);%-ff(:,:,1);
w2=(-rho/dt)*u(:,:,2)+rho*ADV(:,:,2)-mu/2*Laplacian2D(u(:,:,2),N,h);%-ff(:,:,2);
Divw=(w1(xp,:)-w1)/h + (w2(:,yp)-w2)/h;
pp=PoissonSolver2D(Divw,L_hat2);
c3=dt/rho;
r1=u(:,:,1)-dt*ADV(:,:,1)+c2*Laplacian2D(u(:,:,1),N,h)+c3*ff(:,:,1)-c3*(pp-pp(xm,:))/h;
r2=u(:,:,2)-dt*ADV(:,:,2)+c2*Laplacian2D(u(:,:,2),N,h)+c3*ff(:,:,2)-c3*(pp-pp(:,ym))/h;
u1_hat=fft2(r1)./(1-c2*L_hat1);
u2_hat=fft2(r2)./(1-c2*L_hat1);
u(:,:,1)=real(ifft2(u1_hat));
u(:,:,2)=real(ifft2(u2_hat));
ADV=advection2D(u,N,h);

% compute X(n+1)
u0(1)=mean(mean(uu(:,:,1))); u0(2)=mean(mean(uu(:,:,2)));
Cu = Curl2Dedge2node(uu,N,h);
a  = PoissonSolver2D(Cu,L_hat2);
X = X + dt * (interpMAC2Dnew(a,XX,N,h,kID,kID_d,K)+repmat(u0,Ns,1));
areaIBDF(1,2) = polyarea([X(:,1); X(1,1)],[X(:,2); X(1,2)]);
[ampl,ppint,xf,yf] = makeIBcurve(X(:,1)-L/2,X(:,2)-L/2,Ns);
areaIBDF(2,2) = ppint;

% Xtracer = Xtracer + dt * (interpMAC2Dnew(a,XXtracer,N,h,kID,kID_d,K)+repmat(u0,Ntracer,1));
% sk=1;
% areaTracer(1,2) = polyarea([Xtracer(1:sk:end,1); Xtracer(1,1)], [Xtracer(1:sk:end,2);Xtracer(1,2)]);
% [ampl,ppint,xf,yf] = makeIBcurve(Xtracer(1:sk:end,1)-L/2,Xtracer(1:sk:end,2)-L/2,100);
% areaTracer(2,2) = ppint;
% sk=2;
% areaTracer(3,2) = polyarea([Xtracer(1:sk:end,1); Xtracer(1,1)], [Xtracer(1:sk:end,2);Xtracer(1,2)]);
% [ampl,ppint,xf,yf] = makeIBcurve(Xtracer(1:sk:end,1)-L/2,Xtracer(1:sk:end,2)-L/2,100);
% areaTracer(4,2) = ppint;
% sk=20;
% areaTracer(5,2) = polyarea([Xtracer(1:sk:end,1); Xtracer(1,1)], [Xtracer(1:sk:end,2);Xtracer(1,2)]);
% [ampl,ppint,xf,yf] = makeIBcurve(Xtracer(1:sk:end,1)-L/2,Xtracer(1:sk:end,2)-L/2,100);
% areaTracer(6,2) = ppint;

vol(1)=abs(areaIBDF(1,2)-areaIBDF(1,1))/areaIBDF(1,1);
% voltracer(1) = abs(areaTracer(1,2)-areaTracer(1,1))/areaTracer(1,1);

%% main loop
for n = 2 : Nt
    
    uold = u;
    
    % compute X(n+1/2)
    u0(1)=mean(mean(u(:,:,1))); u0(2)=mean(mean(u(:,:,2)));
    Cu = Curl2Dedge2node(u,N,h);
    a  = PoissonSolver2D(Cu,L_hat2);
    XX = X + (dt/2) * (interpMAC2Dnew(a,X,N,h,kID,kID_d,K)+repmat(u0,Ns,1));
    
%      XXtracer = Xtracer + (dt/2) * (interpMAC2Dnew(a,Xtracer,N,h,kID,kID_d,K)+repmat(u0,Ntracer,1));
    
    % compute f(n+1/2)
    FF=Force_HookeSpring(XX,Ns,ds)*ds;
    ff=spreadMAC2Dnew(FF,XX,N,h,kID,kID_d,K,L_hat2);
    
    % solve fluid for u(n+1)
    S=3/2*ADV-1/2*ADVold;
    u = NavierStokes2D_FFT(u,N,h,dt,mu,rho,S,ff,L_hat1,L_hat2);
   
    ADVold = ADV;
    ADV=advection2D(u,N,h);
    
    % compute X(n+1)
    uu = (u+uold)/2;
    uu0(1)=mean(mean(uu(:,:,1))); uu0(2)=mean(mean(uu(:,:,2)));
    Cuu = Curl2Dedge2node(uu,N,h);
    aa  = PoissonSolver2D(Cuu,L_hat2);
    X = X + dt * (interpMAC2Dnew(aa,XX,N,h,kID,kID_d,K)+repmat(uu0,Ns,1));
    
%     Xtracer = Xtracer + dt *  (interpMAC2Dnew(aa,XXtracer,N,h,kID,kID_d,K)+repmat(uu0,Ntracer,1));
    
    
    areaIBDF(1,n+1) = polyarea([X(:,1); X(1,1)],[X(:,2); X(1,2)]);
    [ampl,ppint,xf,yf] = makeIBcurve(X(:,1)-L/2,X(:,2)-L/2,Ns);
    areaIBDF(2,n+1) = ppint;
    vol(n)=abs(areaIBDF(1,n+1)-areaIBDF(1,1))/areaIBDF(1,1);



%     sk=1;
%     areaTracer(1,n+1) = polyarea([Xtracer(1:sk:end,1); Xtracer(1,1)], [Xtracer(1:sk:end,2);Xtracer(1,2)]);
%     [ampl,ppint,xf,yf] = makeIBcurve(Xtracer(1:sk:end,1)-L/2,Xtracer(1:sk:end,2)-L/2,100);
%     areaTracer(2,n+1) = ppint;
%     sk=2;
%     areaTracer(3,n+1) = polyarea([Xtracer(1:sk:end,1); Xtracer(1,1)], [Xtracer(1:sk:end,2);Xtracer(1,2)]);
%     [ampl,ppint,xf,yf] = makeIBcurve(Xtracer(1:sk:end,1)-L/2,Xtracer(1:sk:end,2)-L/2,100);
%     areaTracer(4,n+1) = ppint;
%     sk=20;
%     areaTracer(5,n+1) = polyarea([Xtracer(1:sk:end,1); Xtracer(1,1)], [Xtracer(1:sk:end,2);Xtracer(1,2)]);
%     [ampl,ppint,xf,yf] = makeIBcurve(Xtracer(1:sk:end,1)-L/2,Xtracer(1:sk:end,2)-L/2,100);
%     areaTracer(6,n+1) = ppint;

   
    
    
    %%% caculate fore on the marker
%     Fmarker=Force_HookeSpring(X,Ns,ds)*ds;%*(1+2*0.45*sin(9*(n-0.5)*dt))*4;
%     for k = 1 : Ns
%         [th,rr] = cart2pol(X(k,1)-L/2,X(k,2)-L/2);
%         P = [cos(th), sin(th); -sin(th), cos(th)];
%         Fpol(k,:) = (P * Fmarker(k,:)')'./norm(Fmarker(k,:));
%         TH(k) = th;
%     end

    %% plotting
    if mod(n,floor(Nt/Nf)) == 0
        disp([num2str(n/Nt*100),'%'])
        
        if strcmp(showplot, 'on') 
            
            
            % interp velocity on the marker
            u0(1)=mean(mean(u(:,:,1))); u0(2)=mean(mean(u(:,:,2)));
            Cu = Curl2Dedge2node(u,N,h);
            a  = PoissonSolver2D(Cu,L_hat2);
            Umarker = interpMAC2Dnew(a,X,N,h,kID,kID_d,K)+repmat(u0,Ns,1);
            

            j=n/(Nt/Nf);
            XIBDF(:,:,j+1)=X;
            uIBDF(:,:,:,j+1)=u;
            uu = (u(kp,:,1) + u(:,:,1))/2;
            vv = (u(:,kp,2) + u(:,:,2))/2;
            uMax=max(max(sqrt(uu.^2+vv.^2)));
            UMax=max(sqrt(Umarker(:,1).^2+Umarker(:,2).^2));
            Re=rho*uMax*L/mu;
            %vort = (u(:,:,2)-u(km,:,2) - u(:,:,1) + u(:,km,1) )/h;
                        
            subplot(3,1,[1:2]);
            %contourf(xx,yy,vort ,100); hold on
            %shading flat;  
            %colormap(jet); colorbar; %caxis([-10,10])
            quiver(xx2(sk,sk),yy2(sk,sk),uu(sk,sk),vv(sk,sk),1); hold on;
%             plot(mod(Xtracer(:,1),L),mod(Xtracer(:,2),L),'r-','linewidth',2)
%             plot(mod(Xtracer0(:,1),L),mod(Xtracer0(:,2),L),'g-','linewidth',2)
            plot(mod(X(:,1),L),mod(X(:,2),L),'k-','linewidth',2);            %plot([X(1,1),X(Ns/2,1)],[X(1,2),X(Ns/2,2)],'ro','markerfacecolor','r');
            xlabel(['t = ',num2str(n*dt), ', Re = ',num2str(Re),', N = ',num2str(Nx),...
                    ', # marker per cell = ',num2str(MpC)])
            title(['NS, IB-DF, max|U|= ',num2str(uMax)])
            axis equal
            axis([L/2,.75*L,L/2,.75*L])
            %axis([.2,.5,.5,.8]);
            %axis([.25,.35,.6,.7])
            drawnow
            hold off

%             subplot(6,2,[9,11]);
%             magU= sqrt(uu.^2 + vv.^2);
%             pcolor(xx2,yy2,magU);
%             shading flat;
%             colorbar;
%             axis equal;
%             axis([0,L,0,L])
%             title({'magnitude of velocity: |u|',['max |U|: ',num2str(UMax)]})
%             drawnow

            % plot volume loss
            subplot(3,1,[3]);
%            voltracer(n)=abs(areaTracer(1,n+1)-areaTracer(1,1))/areaTracer(1,1);
            semilogy((1:n)*dt,vol, 'k-','linewidth',1); hold on;
%            semilogy((1:n)*dt,voltracer, 'r-','linewidth',1); 
            ylabel('area loss');
            set(gca,'Ytick',10.^[-14:2:0]);
            axis([0,tend,1e-14,1e-0])
            % legend('IB points', 'Tracers','Location','ne')
            legend('IB points','Location','ne')
            xlabel('t')
            drawnow


            %print('-dpng',['./movies/IBDF1/',num2str(j+10000),'.png'])
        end
    end 
end
                        
    
