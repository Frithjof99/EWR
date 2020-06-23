function uh = fdsolbvp(grd,opr,bcs,pde)
%% CALL: uh = fdsolbvp(grd,opr,bcs,pde);
% INPUT:
%    grd ... STRUCT; Gitter-Struktur.
%    opr ... STRUCT; Operatorinformation.
%    bcs ... STRUCT; Randinformationen.
%    pde ... STRUCT; Koeffizientenfunktionen.
% OUTPUT:
%    uh ... DOUBLE*; Gitterfunktion numerische Loesung.
% DESCRIPTION:
% UH = FDSOLBVP(GRD,OPR,BCS,PDE,VARARGIN) liefert die diskrete Loesung UH des
% Randwertproblems der Eingabe.

% Version 1.0: Willy Doerfler, KIT, 2013.

%% Analysiere Eingabedaten
% Randbedingungen
if nargin < 3 || isempty(bcs)
   bcs = [];
end
% Partielle Differentialgleichung
if nargin < 4 || isempty(pde)
   pde = fdpde(grd,opr,1);
end

%% Reduzierte Matrix aufstellen
[A,b,E,ur] = fdassembvp(grd,opr,bcs,pde);

%% Lösen
uhi = A\b;% Loesung im Inneren
uh  = E*uhi+ur;% Globale Loesung

return
