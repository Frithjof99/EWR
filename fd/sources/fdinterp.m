function vh = fdinterp(grd,opr,func)
%% CALL: vh = fdinterp(grd,opr,func);
% INPUT:
%    grd ... STRUCT; Gitter-Struktur.
%    opr ... STRUCT; Operatorinformation.
%    func ... FUNC; Zu interpolierende Funktion x|->func(x_1,...,x_d).
%                   Cell: {func1,...func_m} -> vh(:,m).
% OUTPUT:
%    vh ... DOUBLE*; vh = [func(x)]_x, x Gitterpunkte.
% DESCRIPTION:
% VH = FDINTERP(GRD,OPR,FUNC) wertet die Funktion FUNC auf dem Rechengitter aus
% und liefert die Gitterfunktion VH zurueck. Ist FUNC vom Typ CELL, so werden
% die Funktionswerte in Spalten von VH ausgegeben.

% Version: Willy Doerfler, KIT, Jun 2020.

% Todo: opr? Effizienz?

%% Analysiere Eingabedaten
dim = grd.dim;

%% Vektor belegen
if dim==1
   ind = grd.G>0;% Unindizierte Punkte ausfiltern
   G = grd.G(ind);
   x = grd.x(ind);
   vh = zeros(length(ind),1);
   vh(G) = func(x);
elseif dim==2
   if ~iscell(func)% Funktionen, skalar- oder vektorwertig
      vt = func(grd.x(1),grd.y(1));% Testausgabe
      vh = zeros(max(max(grd.G)),size(vt,2));
      if numel(grd.G)==nnz(grd.G)% Keine Nullen: vermeidet 'if G>0'
         for i=1:length(grd.x)
            for j=1:length(grd.y)
               vh(grd.G(i,j),:) = func(grd.x(i),grd.y(j));
            end
         end
      else
         for i=1:length(grd.x)
            for j=1:length(grd.y)
               if grd.G(i,j)>0
                  vh(grd.G(i,j)) = func(grd.x(i),grd.y(j));
               end
            end
         end
      end
   else% Zelle skalarwertiger Funktionen
      lf = length(func);
      vh = zeros(max(max(grd.G)),lf);
      if numel(grd.G)==nnz(grd.G)% Keine Nullen: vermeidet 'if G>0'
         for kk=1:lf
            for i=1:length(grd.x)
               for j=1:length(grd.y)
                  vh(grd.G(i,j),kk) = func{kk}(grd.x(i),grd.y(j));
               end
            end
         end
      else
         for kk=1:lf
            for i=1:length(grd.x)
               for j=1:length(grd.y)
                  if grd.G(i,j)>0
                     vh(grd.G(i,j),kk) = func{kk}(grd.x(i),grd.y(j));
                  end
               end
            end
         end
      end      
   end
end

return
