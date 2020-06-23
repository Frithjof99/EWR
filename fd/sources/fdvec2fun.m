function y = fdvec2fun(y,grd)
%% Call: y = fdvec2fun(y,grd);
% Input:
%    y ... sparse DBLE(*,1); Nodal vector as in grd.
% Output:
%    y ... sparse DBLE(*,1); Solution vector lexicographic.
% Description:
%%
y(grd.G) = y;
return
%%
