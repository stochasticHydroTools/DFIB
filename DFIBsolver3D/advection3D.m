function adv = advection3D(u,N,h)

xp=[2:N(1),1];
xm=[N(1),1:N(1)-1];

yp=[2:N(2),1];
ym=[N(2),1:N(2)-1];

zp=[2:N(3),1];
zm=[N(3),1:N(3)-1];

u1 = u(:,:,:,1);
u2 = u(:,:,:,2);
u3 = u(:,:,:,3);


% 1st component of advection
Dxu1 = (u1(xp,:,:) - u1(xm,:,:)) / (2*h);
Dyu1 = (u1(:,yp,:) - u1(:,ym,:)) / (2*h);
Dzu1 = (u1(:,:,zp) - u1(:,:,zm)) / (2*h);
u2i  = (u2+u2(xm,:,:))/2; u2i = (u2i+u2i(:,yp,:))/2;
u3i  = (u3+u3(xm,:,:))/2; u3i = (u3i+u3i(:,:,zp))/2;
u1u1 = u1.*u1;  Dxu1u1 = (u1u1(xp,:,:)-u1u1(xm,:,:))/(2*h);
u1u2 = u1.*u2i; Dyu1u2 = (u1u2(:,yp,:)-u1u2(:,ym,:))/(2*h);
u1u3 = u1.*u3i; Dzu1u3 = (u1u3(:,:,zp)-u1u3(:,:,zm))/(2*h);
adv(:,:,:,1) = 0.5 * (u1.*Dxu1+u2i.*Dyu1+u3i.*Dzu1) + ...
               0.5 * (Dxu1u1 + Dyu1u2 + Dzu1u3);
           
% 2nd component of advection
u1i = (u1+u1(:,ym,:))/2; u1i = (u1i+u1i(xp,:,:))/2;
u3i = (u3+u3(:,ym,:))/2; u3i = (u3i+u3i(:,:,zp))/2;
Dxu2= (u2(xp,:,:)-u2(xm,:,:))/(2*h);
Dyu2= (u2(:,yp,:)-u2(:,ym,:))/(2*h);
Dzu2= (u2(:,:,zp)-u2(:,:,zm))/(2*h); 
u1u2= u1i.*u2; Dxu1u2 = (u1u2(xp,:,:)-u1u2(xm,:,:))/(2*h);
u2u2= u2.*u2 ; Dyu2u2 = (u2u2(:,yp,:)-u2u2(:,ym,:))/(2*h);
u2u3= u2.*u3i; Dzu2u3 = (u2u3(:,:,zp)-u2u3(:,:,zm))/(2*h);
adv(:,:,:,2) = 0.5 * (u1.*Dxu2+u2i.*Dyu2+u3i.*Dzu2) + ...
               0.5 * (Dxu1u2 + Dyu2u2 + Dzu2u3);
           
% 3rd compoent of advection
u1i = (u1+u1(:,:,zm))/2; u1i = (u1i+u1i(xp,:,:))/2;
u2i = (u2+u2(:,:,zm))/2; u2i = (u2i+u2i(:,yp,:))/2;
Dxu3= (u3(xp,:,:)-u3(xm,:,:))/(2*h);
Dyu3= (u3(:,yp,:)-u3(:,ym,:))/(2*h);
Dzu3= (u3(:,:,zp)-u3(:,:,zm))/(2*h);
u1u3= u1i.*u3; Dxu1u3 = (u1u3(xp,:,:)-u1u3(xm,:,:))/(2*h);
u2u3= u2i.*u3; Dyu2u3 = (u2u3(:,yp,:)-u2u3(:,ym,:))/(2*h);
u3u3= u3.*u3 ; Dzu3u3 = (u3u3(:,:,zp)-u3u3(:,:,zm))/(2*h);
adv(:,:,:,3) = 0.5 * (u1.*Dxu3 + u2i.*Dyu3 + u3i.*Dzu3) +...
               0.5 * (Dxu1u3 + Dyu2u3 + Dzu3u3);