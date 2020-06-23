function opr = fdopr(grd,s1,stc,s2,order1)
%% CALL: opr = fdopr(grd,s1,stc,s2,order1);
% INPUT:
%    grd ... STRUCT; Gitter-Struktur.
%    stc ... INT; Stencilsize in points [DEF dim=1: 3, dim=2: 5].
% OUTPUT:
%    opr ... STRUCT; Operatorinformation.
% DESCRIPTION:
% FDOPR Information ueber den Differentialoperator.
%    stencil ... DIM=1: 3pt; DIM=2: 5pt, 9pt; DIM=3: 7pt
%    order1  ... 1. Ordnung Diff.-op: Typ '', 'central' oder 'upwind'

% Version 2.0: Willy Doerfler, KIT, 2020.

%% Initialisierungen
dim = grd.dim;

%% Analysiere Input
if nargin < 5
   if dim==1, stc = 3; end
   if dim==2, stc = 5; end
   if dim==3, stc = 7; end
   order1 = 'central';
end

%% Setzen
opr.stencil = stc;
opr.order1  = order1;
%
if dim==3% Spezialfall
   opr = fdopr3(grd,s1,stc,s2,order1);
end

return
