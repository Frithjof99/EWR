function [X,Y,Z] = fdmeshgrid(grd,fh)
%% CALL: [X,Y,Z] = fdmeshgrid(grd,fh);
% INPUT:
%    grd ... STRUCT; Gitter-Struktur.
%    fh ... DOUBLE*/FUNC; Gitterfunktion [DEF 0].
% OUTPUT:
%    X,Y,Z ... DOUBLE**; Meshgrid Daten.
% DESCRIPTION:
% FDMESHGRID Bereitet die Visualisierung von Gitterfunktionen vor.
% [X,Y,Z] = FDMESHGRID(GRD,FH) gibt Matrizen X, Y und Z zurück, welche als
% Eingabeparameter für die Funktionen SURF, MESH oder CONTOUR verwendet werden
% können. Das Finite-Differenzen-Gitter ist durch die Struktur GRD mit den
% Komponenten G, X, ... definiert (siehe FDGRID). Die Werte der Gitterfunktion
% sind durch eine Funktion x|->fh(x_1,...,x_d) oder den Vektor FH definiert.
%   FH Funktionspointer: Auswertung fuer G>0,
%   FH Vektor(N):  Auswertung fuer G>0.
%   FH Vektor(NI): Auswertung nur in inneren Knoten.

% Version 0.0: Markus Richter, KIT, 2008.
% Willy Doerfler, KIT, Jun 2020.

%% Initialisierungen
[m,n] = size(grd.G);
X = nan(m,n);% Macht den Rest unsichtbar
Y = nan(m,n);

%% Erzeuge die X-, Y-, und Z-Matrizen für die Funktionen surf, mesh, etc.  
if isnumeric(fh)% fh ist Vektor
   if size(fh,1)==grd.ni% Ggf. Rand ausblenden
      i = grd.bverts(1,:);% Randknoten
      j = grd.bverts(2,:);
      grd.G(i,j) = 0;
   end
   kk = size(fh,2);% Spalten
   if kk==1% Eine Spalte
      Z = nan(m,n);
      for i=1:m
         for j=1:n
            if grd.G(i,j)>0
               X(i,j) = grd.x(i);
               Y(i,j) = grd.y(j);
               Z(i,j) = fh(grd.G(i,j));
            end
         end
      end
   else% Mehr Spalten
      Z = nan(m,n,kk);
      for i=1:m
         for j=1:n
            if grd.G(i,j)>0
               X(i,j) = grd.x(i);
               Y(i,j) = grd.y(j);
               for ll=1:kk
                  Z(i,j,ll) = fh(grd.G(i,j),ll);
               end
            end
         end
      end
   end
else% fh ist Funktionspointer
   Z = nan(m,n);
   for i=1:m
      for j=1:n
         if grd.G(i,j)>0
            X(i,j) = grd.x(i);
            Y(i,j) = grd.y(j);
            Z(i,j) = fh(grd.x(i),grd.y(j));
         end
      end
   end
end

return
