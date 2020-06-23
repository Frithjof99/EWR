function [bcs,grd] = fdbcs(grd,essbpts,essbfuns,natbpts,natbfuns,perbpts,perbfuns)
%% CALL: [bcs,grd] = fdbcs(grd,essbpts,essbfuns,natbpts,natbfuns,perbpts,perbfuns);
% INPUT:
%    grd ... STRUCT; Gitter-Struktur.
%    essbpts ... FUNC; True fuer Dirichletrandpunkte [Def].
%    essbfuns ... FUNC; Dirichletrandwerte-Funktion [Def 0].
%    natbpts ... FUNC; True fuer Neumannrandpunkte.
%    natbfuns ... FUNC; Neumannrandwerte-Funktion.
%    perbpts ... FUNC; True fuer Punkte des periodischen Randes (???).
%    perbfuns ... FUNC; uP (s.u.).
% OUTPUT:
%    bcs ... STRUCT; Randinformationen.
% DESCRIPTION:
% FDBCS Randpunkte und Randwertfunktionen.
%    essb ... INT; Flag fuer Dirichletrand.
%    natb ... INT; Flag fuer Neumannrand.
%    perb ... INT; Flag fuer periodischen Rand.
%    essbfuns ... Essentielle Randbedingung (Dirichlet): u = uD.
%    natbfuns ... Natuerliche Randbedingung (Neumann): a d_n u = uN.
%    perbfuns ... Periodische Randbedingung: u(x) = uP(u(xP)).

% ToDo: Periodische Randbedingungen

% Version 1.0: Willy Doerfler, KIT, 2013.
% Version 3.0: Willy Doerfler, KIT, 2017.

%% Analysiere Eingabe
if nargin<2
   essbpts = 1; essbfuns = 0; natbpts = []; natbfuns = [];
   perbpts = []; perbfuns = [];
elseif nargin<3
   essbfuns = []; natbpts = []; natbfuns = []; perbpts = []; perbfuns = [];
elseif nargin<4
   natbpts = []; natbfuns = []; perbpts = []; perbfuns = [];
elseif nargin<5
   natbfuns = []; perbpts = []; perbfuns = [];
elseif nargin<6
   perbpts = []; perbfuns = [];
elseif nargin<7
   perbfuns = [];
end

%% Initialisieren
essb = 1; natb = 1; perb = 1;
if isempty(essbpts), essb = 0; end
if isempty(natbpts), natb = 0; end
if isempty(perbpts), perb = 0; end

%% Besetzen
bcs.essb = essb;
bcs.natb = natb;
bcs.perb = perb;
bcs.essbpts = essbpts;
bcs.natbpts = natbpts;
bcs.perbpts = perbpts;
bcs.essbfuns = essbfuns;
bcs.natbfuns = natbfuns;
bcs.perbfuns = perbfuns;

%% Randinformation in grd bei Bedarf setzen
if nargout==2% Bedarf
   if essb% Dirichlet Bedingung aktiv
      if grd.dim==1
         grd.bverts(2,:) = ones(1,size(grd.bverts,2));% Noch global
      end
      if grd.dim==2
         grd.bverts(3,:) = ones(1,size(grd.bverts,2));
      end
   end
end

return
