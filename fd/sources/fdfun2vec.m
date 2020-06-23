function y = fdfun2vec(y,grd)
%% Call: y = fdfun2vec(y,grd);
% Input:
%    y ... sparse DBLE(*,1); Solution vector lexicographic.
% Output:
%    y ... sparse DBLE(*,1); Nodal vector as in grd.
% Description:
%%
y = y(grd.G);
return
%%
