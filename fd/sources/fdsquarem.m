function grd = fdsquarem(n,geom)
%% CALL: grd = fdsquarem(n,geom);
% INPUT:
%    n ... INT; Gitterpunkte pro dimension [DEF n=10].
%    geom ... DOUBLE; Geometrie des Rechengebietes [DEF [zeros(dim) 2]].
% OUTPUT:
%    grd ... STRUCT; Gitter-Struktur.
% DESCRIPTION:
% FDSQUAREM Finite-Differenzen-Gitter auf einem Quadrat.
% Beschreibung GEOM:
%     [c size] - 'c' center in IR^d, 'size' Durchmesser

% Version 1.0: Willy Doerfler, KIT, 2012.

%% Analysiere Input
switch nargin
case 0
   n = 10;
   geom = 2;
case 1
   geom = 2;
end

%% Initialisierungen
if length(geom)>1
   r = geom;
else
   r = geom*[0,0,1];
end

%% Gitter holen
% Transf. [-1,1]^2->[r(1),r(2)]x[r(3),r(4)]
grd = fdgrid('S',n);
grd = fdtransform(grd,[r(1);r(2)],r(3)/2*eye(2));

return
