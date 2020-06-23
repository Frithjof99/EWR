function bcs = fddirichlet(grd,uD)
%% CALL: bcs = fddirichlet(grd,uD);
% INPUT:
%    grd ... STRUCT; Gitter-Struktur.
%    uD ... FUNC; Randwertefunktion fuer Dirichletrandbedingung.
% OUTPUT:
%    bcs ... STRUCT; Randinformationen.
% DESCRIPTION:
% FDDIRICHLET Randwertefunktion fuer Dirichletrandbedingung auf dem ganzen Rand.
%    Setzt Randbedingung: u(x) = uD(x). [DEF 0]
% Verwende [] fuer nichtaktive Funktion.

% Version 1.0: Willy Doerfler, KIT, 2013.
% Version 3.0: Willy Doerfler, KIT, 2017.

switch nargin
case 1
   switch grd.dim
      case 1, bcs = fdbcs(grd,1,@(xx)0);
      case 2, bcs = fdbcs(grd,1,@(xx,yy)0);
      otherwise, error('*** Error *** dim\not\in{1,2} in fddirichlet');
   end
case 2
   bcs = fdbcs(grd,1,uD);
otherwise
   error('*** error *** Check syntax for fddirichlet');
end

return
