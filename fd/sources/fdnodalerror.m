function err = fdnodalerror(grd,vh,exakt)
%% CALL: err = fdnodalerror(grd,vh,exakt);
% INPUT:
%    grd ... STRUCT; Gitter-Struktur.
%    vh ... DOUBLE*; Gitterfunktion [DEF 0].
%    exakt ... FUNC; Exakte Loesung als Funktion von x_1,...,x_d.
% OUTPUT:
%    err ... DOUBLE; max_{x in G}|vh(x)-exakt(x)|.
% DESCRIPTION:
% FDNODALERROR berechnet den nodalen Fehler.
% ERR = FDNODALERROR(GRD,VH,EXAKT) wertet die maximale Differenz zwischen dem
% Vektor VH auf dem Gitter GRD  und der Funktion EXAKT in den Gitterpunkten aus
% GRD aus.
% Ist LENGTH(VH)=GRD.NN, geschieht dies auf allen Punken, fuer LENGTH(VH)=GRD.NI
% nur auf den inneren Punkten.

% Version 1.0: Willy Doerfler, KIT, 2013.
% Version 2.0: Willy Doerfler, KIT, 2020.

%% Analysiere Eingabedaten
if length(vh)==grd.ni% Ggf. Rand ausblenden
   ir = grd.bverts(1,:);% Randknoten
   jr = grd.bverts(2,:);
   grd.G(ir,jr) = 0;
end

%% Maximalfehler
[m,n] = size(grd.G);
err = 0;
if grd.dim==1
   for i=1:n
      gi = grd.G(i);
      if gi>0
         err = max(err,abs(exakt(grd.x(i))-vh(gi)));
      end
   end
elseif grd.dim==2
   for i=1:m
      for j=1:n
         gij = grd.G(i,j);
         if gij>0
            err = max(err,abs(exakt(grd.x(i),grd.y(j))-vh(gij)));
         end
      end
   end
elseif grd.dim==3
   err = fdnodalerror3(grd,vh,exakt,m);
end

return
