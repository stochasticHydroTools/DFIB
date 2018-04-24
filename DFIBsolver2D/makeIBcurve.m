%-------------------------------------------------------------------------%
% makeIBcurve.m
% Yuanxun Bao
% Mar 3, 2014
% 
% output a parametric curve of give IB points
% x  - x-coord of IB points (curve c)
% y  - y-coord of IB points
% xf - x-coord of IB curve
% yf - y-coord of IB curve
%-------------------------------------------------------------------------%

function [ampl,ppint,xf,yf] = makeIBcurve(x,y,Nf)

% convert to polar
[th,r]=cart2pol(x,y);
% sort th in increasing order
[th,index]=sort(th); r=r(index);

% interpolate with cubic spline with periodic BCs
% pp is a parametric curve
pp = csape([th;th(1)+2*pi],[r;r(1)],'perodic');

% eval on equally spaced grid [-pi pi]
th_f = linspace(-pi,pi,Nf);
th_f = (th_f < th(1)).*(th_f+2*pi) + ...
       (th_f>=th(1) & th_f<=th(1)+2*pi) .* th_f + ...
       (th_f > th(1)+2*pi) .* (th_f-2*pi);
r_f  = ppval(pp, th_f);

r_f_hat = fft(r_f(1:end-1)-1);
ampl = r_f_hat(3)*2; % the 3rd index is associated with cos(2*theta) mode
ampl = real(ampl)/(Nf-1);

% convert back to cartesian
[xf,yf] = pol2cart(th_f,r_f);
% xf = r_f;
% yf = th_f;

% find the area enclosed by the curve: (1/2)*int(r(th)^2, th = a..b )
ppint = 0;
for k = 1 : pp.pieces
    h  = pp.breaks(k+1)-pp.breaks(k);
    a  = pp.coefs(k,:);
    ppint = ppint + 1/7*a(1)^2*h^7 + 1/3*a(1)*a(2)*h^6+1/5*(2*a(1)*a(3)+a(2)^2)*h^5 ...
          + 1/4*(2*a(1)*a(4)+2*a(2)*a(3))*h^4 + 1/3*(2*a(2)*a(4)+a(3)^2)*h^3 ...
          + a(3)*a(4)*h^2+a(4)^2*h; 
end
ppint = ppint/2;