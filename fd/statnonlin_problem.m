%% Zusammenfassung: Stationaeres nichtlineares Randwertproblem
% Version: Willy Doerfler, KIT, Jun 2020.

%% Gitter erzeugen
grd = fdgrid(domain,NperDim,geom);

%% Gebiet darstellen
figure(1); clf
   fdplot(grd);
   title(['Rechengebiet (' num2str(grd.nn) ')']);
drawnow

%% Konfiguration des FD-Operator
opr = fdopr(grd,'stencil',star);

%% PDE definieren (Koffizientenfunktionen)
% Allgemeine Form:
%    - div(a grad u) + b.grad u + cu = nlf(u)   in Omega,
%                                  u = dir      auf rand(Omega).
% Zudem wird nlf' = nlfp fuer das Newtonverfahren benoetigt.
pde = fdpde(grd,opr,fun_a,fun_b,fun_c,fun_f);
bcs = fddirichlet(grd,fun_dir);

%% Reduzierte Matrix aufstellen
[A,b,E,ur] = fdassembvp(grd,opr,bcs,pde);

%% Loesen

% Initialisierung
u0 = ur;% Startwert u0
v0 = E'*(u0-ur);% u0 zu v0 reduzieren
NlFunc  = @(vv) A*vv-b-E'*fun_nlf(E*vv+ur);
NlFuncp = @(vv) A-E'*spdiags(fun_nlfp(E*vv+ur),0:0,grd.nn,grd.nn)*E;

% Newton
nopt = newton(); nopt.Display = 3; nopt.StepCtrl = 1;% Newton-Optionen
[uhi,fu] = newton(NlFunc,NlFuncp,v0,nopt);
uh = E*uhi+ur;% Globale Loesung
fprintf(' Newtonresidual: %7.4e\n',norm(fu,2));

%% Loesung graphisch darstellen
figure(2); clf
   fdplot(grd,uh,gopt);
   t = title(['Numerische Loesung (lambda = ' num2str(lambda) ')']);
   t.FontSize = 25;
drawnow

%% Ggf Fehler berechnen
if ~isempty(fun_uex)
   errli = fdnodalerror(grd,uh,fun_uex);
   errl2 = fdl2error(grd,uh,fun_uex);
   fprintf(' Anzahl Freiheitsgrade: %d\n',grd.ni);
   fprintf(' Nodaler Fehler: %7.4e\n',errli);
   fprintf(' l2-Fehler     : %7.4e\n',errl2);
end

%% ENDE
