function f = Laplacian3D(phi,N,h)


xp = [2:N(1),1];
xm = [N(1),1:N(1)-1];

yp = [2:N(2),1];
ym = [N(2),1:N(2)-1];

zp = [2:N(3),1];
zm = [N(3),1:N(3)-1];

f = (phi(xp,: ,: ) + phi(xm,: ,: ) + ...
     phi(: ,yp,: ) + phi(: ,ym,: ) + ...
     phi(: ,: ,zp) + phi(: ,: ,zm) - 6*phi) / (h*h);