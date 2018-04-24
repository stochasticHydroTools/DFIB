function v = fluidSolverMAC(w,nu,rho,dt,N,h,L_hat)

kp = [2:N,1];
km = [N,1:N-1];

Dw = ( w(kp,:,1) - w(:,:,1) + w(:,kp,2) - w(:,:,2) ) / h;
Dw_hat = fft2(Dw);

q_hat = Dw_hat ./ L_hat;
q = real(ifft2(q_hat));

rhs1 = w(:,:,1) - (q - q(km,:))/h;
rhs2 = w(:,:,2) - (q - q(:,km))/h;

c = (nu*dt)/(2*rho);

L_hat(1,1) = 0; % don't need L_hat(1,1)=1.
v1_hat = fft2(rhs1) ./ (1-c*L_hat); 
v2_hat = fft2(rhs2) ./ (1-c*L_hat); 

v(:,:,1) = real(ifft2(v1_hat));
v(:,:,2) = real(ifft2(v2_hat));