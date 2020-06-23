function [E,ur] = fdassembcs(grd,opr,bcs)
%% CALL: [E,ur] = fdassembcs(grd,opr,bcs);
% INPUT:
%    grd ... STRUCT; Gitter-Struktur.
%    opr ... STRUCT; Operatorinformation (Def []).
%    bcs ... STRUCT; Randinformationen.
% OUTPUT:
%    E ... DOUBLE**; E ist die Fortsetzung des reduzierten Vektors auf den
%                    globalen Vektor.
%    ur ... DOUBLE*; Vektor mit Dirichlet-Randwerten.
% DESCRIPTION:
% FDASSEMBCS Dirichlet-Randwerte in einem Vektor sammeln und Erweiterungsmatrix.
% [E,UR] = FDASSEMDIR(GRD,DIRFUN) gibt den globalen Vektor mit den
% Dirichlet-Randwerten und mit Werten 0 im Inneren zurueck. Der globale Vektor
% ist der Vektor aller Knotenwerte, im reduzierten Vektor sind die Knoten mit
% essentiellen Randbedingungen entfernt.

% ToDo: Periodische Randbedingungen

% Version 1.0: Willy Doerfler, KIT, 2020.

%% Initialisierungen
nn = grd.nn;
ni = grd.ni;
nb = nn-ni;
ur = zeros(nn,1);
jb = zeros(nb,1);

%% Randwerte belegen
if grd.dim==1
   if ~isempty(bcs.essb)% Dirichletrandpunkte vorhanden
      dfun = bcs.essbfuns;
      % Faelle: func=[] oder func=const zu func
      if isempty(dfun),       dfun = @(xx) 0;
      elseif isnumeric(dfun), dfun = @(xx) dfun; end
      for ir=1:nb
         j = grd.bverts(1,ir);
         ur(j) = dfun(grd.x(j));
         jb(ir) = j;
      end
      % Einschraenkungsmatrix
      E = speye(nn);
      E(:,jb) = [];% Loesche Randknoten
   end
   if bcs.perb~=0% Periodische Randbedingung (Torus)
      % Einschraenkungsmatrix
      E = speye(nn);
      % '1' ist Master, 'nn' ist Slave
      E(nn,1) = 1;% Slave-Knoten eintragen
      E(:,nn) = [];% Loesche Slave-Knoten
   end
elseif grd.dim==2
   if ~isempty(bcs.essb)% Dirichletrandpunkte vorhanden
      dfun = bcs.essbfuns;
      % Faelle: func=[] oder func=const zu func
      if isempty(dfun),       dfun = @(xx,yy) 0;
      elseif isnumeric(dfun), dfun = @(xx,yy) dfun; end
      for ir=1:nb
         i = grd.bverts(1,ir);
         j = grd.bverts(2,ir);
         k = grd.G(i,j);
         ur(k) = dfun(grd.x(i),grd.y(j));
         jb(ir) = k;
      end
      % Einschraenkungsmatrix
      E = speye(nn);
      E(:,jb) = [];% Loesche Randknoten
   end
   if ~isempty(bcs.perb)% Periodische Randbedingung (Torus)
      % Einschraenkungsmatrix
      E = speye(nn);
      for ir=1:nb
         rnr = grd.bverts(3,ir);
         if rnr>0
            slave = grd.G(grd.bverts(1,ir) ,grd.bverts(2,ir));
            %maste = grd.G(grd.bverts(1,rnr),grd.bverts(2,rnr));
            jb(ir) = slave;
            %E(slave,maste) = 1;% Identifikation
         end
      end
      E(:,jb(jb>0)) = [];% Loesche Slave-Knoten
   end
elseif grd.dim==3
   [E,ur] = fdassembcs3(grd,opr,bcs);
end

return
