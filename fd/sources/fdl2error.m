function err = fdl2error(grd,vh,exakt)
%% CALL: err = fdl2error(grd,vh,exakt);
% INPUT:
%    grd ... STRUCT; Gitter-Struktur.
%    vh ... DOUBLE*; Gitterfunktion [DEF 0].
%    exakt ... FUNC; Exakte Loesung als Funktion von x_1,..,x_d.
% OUTPUT:
%    err ... DOUBLE; sqrt( sum_{x in G}|vh(x)-exakt(x)|^2 h_1...h_d ).
% DESCRIPTION:
% FDL2ERROR berechnet den l2-Fehler.
% ERR = FDL2ERROR(GRD,VH,EXAKT) wertet die Summe der quadratischen Differenzen
% zwischen dem Vektor VH auf dem Gitter GRD und der Funktion EXAKT in den
% Gitterpunkten aus GRD aus.
% Ist LENGTH(VH)=GRD.NN geschieht dies auf allen Punken, fuer LENGTH(VH)=GRD.NI
% nur auf den inneren Punkten.

% Version 1.0: Willy Doerfler, KIT, 2013.

%% Analysiere Eingabedaten
if length(vh)==grd.ni% Ggf. Rand ausblenden
   ir = grd.bverts(1,:);% Randknoten
   jr = grd.bverts(2,:);
   grd.G(ir,jr) = 0;
end

%% l2-Fehler
[m,n] = size(grd.G);
err = 0;
if grd.dim==1
   for i=1:n
      gi = grd.G(i);
      if gi>0
         err = err+(exakt(grd.x(i))-vh(gi))^2;
      end
   end
   h = grd.x(2)-grd.x(1);
   err = realsqrt(h*err);
elseif grd.dim==2
   for i=1:m
      for j=1:n
         gij = grd.G(i,j);
         if gij>0
            err = err+(exakt(grd.x(i),grd.y(j))-vh(gij))^2;
         end
      end
   end
   hx = grd.x(2)-grd.x(1);
   hy = grd.y(2)-grd.y(1);
   err = realsqrt(hx*hy*err);
elseif grd.dim==3
   err = fdl2error3(grd,vh,exakt);
end

return
