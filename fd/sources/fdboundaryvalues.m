function ur = fdboundaryvalues(grd,dfun)
%% CALL: ur = fdboundaryvalues(grd,dfun);
% INPUT:
%    grd ... STRUCT; Gitter-Struktur.
%    dfun ... STRUCT; Randwerte.
% OUTPUT:
%    ur ... DOUBLE*; Vektor mit Dirichlet-Randwerten, 0 sonst.
% DESCRIPTION:
% FDBOUNDARYVALUES Dirichlet-Randwerte in einem Vektor sammeln, globaler Vektor
% mit den Dirichlet-Randwerten und mit Werten 0 im Inneren.

% Version 1.0: Willy Doerfler, KIT, Jun 2020.

%% Initialisierungen
nn = grd.nn;
ni = grd.ni;
nb = nn-ni;

%% Randwerte belegen
if grd.dim==1
   % Faelle: dfun=[] oder dfun=const zu dfun
   if isempty(dfun),       dfun = @(xx) 0;
   elseif isnumeric(dfun), dfun = @(xx) dfun; end
   ur = zeros(nn,1);
   for ir=1:nb
      j = grd.bverts(1,ir);
      ur(j) = dfun(grd.x(j));
   end
elseif grd.dim==2
   % Faelle: dfun=[] oder dfun=const zu dfun
   if isempty(dfun),       dfun = @(xx,yy) 0;
   elseif isnumeric(dfun), dfun = @(xx,yy) dfun; end
%    for ir=1:nb
%       i = grd.bverts(1,ir);
%       j = grd.bverts(2,ir);
%       k = grd.G(i,j);
%       ur(k) = dfun(grd.x(i),grd.y(j));
%    end
   % Ziel: Vektorisierter Funktionsaufruf
   k  = zeros(nb,1);
   xr = zeros(nb,1);
   yr = zeros(nb,1);
   for ir=1:nb
      i = grd.bverts(1,ir);
      j = grd.bverts(2,ir);
      k(ir)  = grd.G(i,j);
      xr(ir) = grd.x(i);
      yr(ir) = grd.y(j);
   end
   uu = dfun(xr(1),yr(1));% Test output
   su = size(uu,2);
   if su==1
      ur = zeros(nn,1);
      ur(k(:)) = dfun(xr,yr);
   else
      ur = zeros(nn,su);
      ur(k(:),:) = dfun(xr,yr);
   end
elseif grd.dim==3
   [~,ur] = fdassembcs3(grd,opr,bcs);
end

return
