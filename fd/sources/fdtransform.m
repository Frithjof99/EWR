function grd = fdtransform(grd,bv,Bm)
%% CALL: grd = fdtransform(grd,bv,Bm);
% INPUT:
%    grd ... STRUCT; Gitter-Struktur (Raumdimension dim).
%    bv ... DOUBLE(dim); Verschiebung des Gebietes [DEF zero(dim,1)].
%    Bm ... DOUBLE(dim,dim); Drehung des Gebietes [DEF eye(dim)].
% OUTPUT:
%    grd ... STRUCT; Gitter-Struktur mit transformierten Gitterpunkten.
% DESCRIPTION:
% FDTRANSFORM Lineare Transformation des Rechengebietes.
% GRD = FDTRANSFORM(GRD,BV,BM) transformiert die Koordinaten der 
% Gitterstruktur GRD mittles einer Bewegung X |-> BM*X+BV. In diesem
% Zusammenhang ist nur ein diagonales BM sinnvoll.

% Version 1.0: Willy Doerfler, KIT, 2013.

switch grd.dim
%% DIM=1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 1   
   if nargin < 2
      return
   end
   if nargin < 3
      Bm = 1;
   end
   x = Bm*grd.x + bv;
   grd.x = x;

%% DIM=2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 2

   %% Analysiere Eingabedaten
   if nargin < 2
      return
   end
   if nargin < 3
      Bm = eye(2);
   end
   Bm = diag(diag(Bm));% Nur Diaganalanteil nehmen
  
   %% Transformieren
   x = Bm(1,1)*grd.x + Bm(1,2)*grd.y + bv(1);
   y = Bm(2,2)*grd.y + Bm(2,1)*grd.x + bv(2);
   grd.x = x; 
   grd.y = y;

otherwise
   error(' *** Error *** grd.dim out of range');
end

return
