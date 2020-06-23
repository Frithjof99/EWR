function bcs = fdperiodic(grd,uP)
%% CALL: bcs = fdperiodic(grd,uP);
% INPUT:
%    grd ... STRUCT; Gitter-Struktur.
%    uP ... FUNC; Randwertefunktion fuer periodische Randbedingung (s.u.).
% OUTPUT:
%    bcs ... STRUCT; Randinformationen.
% DESCRIPTION:
% FDPERIODIC Randwertefunktion fuer periodische Randbedingung. Setzt periodische
% Randbedingung: u(x) = uP(u(xP)) [DEF id], xP ist der identifizierte Randpunkt.
% Verwende [] fuer nichtaktive Funktion.

% Version 1.0: Willy Doerfler, KIT, 2013.

switch nargin
case 2
   bcs = fdbcs(grd,[],[],[],[],1,@(uu)uu);
case 3
   bcs = fdbcs(grd,[],[],[],[],1,uP);
otherwise
   error('*** Error *** Check syntax for fdperiodic');
end

return

