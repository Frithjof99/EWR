function ud = fun_u_sing(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call: ud = fun_u_sing(x);
% Input:
%    x ... Nx2; Vector of points
% Output:
%    ud ... Nx1; Vector of function values
% Description:
% Implements the harmonic function
%    u(x) = |x|^{2/3}*cos(2/3\phi(x)-pi/6)
% on a corner domain with vertex=0 and angle in -1/2*pi,...,pi.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r = x(:,1).*x(:,1)+x(:,2).*x(:,2);% r is in fact r^2 !
argx = atan2(x(:,2),x(:,1));
ud = r.^(1/3).*cos((2/3)*argx-pi/6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%