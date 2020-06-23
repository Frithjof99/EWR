function pde = fdpde( ~ , ~ ,stifcoefs,trnscoefs,masscoefs,loadfuns)
%% CALL: pde = fdpde(grd,opr,stifcoefs,trnscoefs,masscoefs,loadfuns);
% INPUT:
%    grd ... STRUCT; Gitter-Struktur.
%    opr ... STRUCT; Operatorinformation.
%    stifcoefs, etc ... FUNC/DOUBLE*;
%                       Beschreibung Koeffizientenfunktionen a,b,c (s.u.).
% OUTPUT:
%    pde ... STRUCT; Koeffizientenfunktionen.
% DESCRIPTION:
% FDPDE Packt die Koeffizientenfunktionen der linearen partiellen 
% Differentialgleichung
%    - a \Delta u + b . \grad u + c u = f .
% in eine Struktur. Die Koeffizientenfunktionen seien von der Form (zB)
% a:[x_1,...,x_d]->IR. a=[] und a=konst sind zulaessig. 

% Version 0.0: Markus Richter, KIT, 2008.
% Version 1.0: Willy Doerfler, KIT, 2013.

%% Analysiere Eingabedaten
% Steifigkeitskoeffizient 'a'
if nargin < 3
   stifcoefs = [];
end
% Transportkoeffizient 'b'
if nargin < 4
   trnscoefs = [];
end
% Reaktionskoeffizient 'c'
if nargin < 5
   masscoefs = [];
end
% Quelldichte 'f'
if nargin < 6
   loadfuns = []; 
end

%% Fuelle Struktur
pde.stifcoefs = stifcoefs;
pde.trnscoefs = trnscoefs;
pde.masscoefs = masscoefs;
pde.loadfuns  = loadfuns;

return
