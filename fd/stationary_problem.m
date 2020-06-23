%% Zusammenfassung: Stationaeres Randwertproblem
% Version: Willy Doerfler, KIT, Jun 2020.

%% Gitter erzeugen
grd = fdgrid(domain,NperDim,geom);

%% Gebiet darstellen
figure(1); clf
   fdplot(grd);
   title(['Rechengebiet (' num2str(grd.nn) ')']);
drawnow

%% Konfiguration des FD-Operator
opr = fdopr(grd,'stencil',star,'order1',order1);

%% PDE definieren (Koffizientenfunktionen)
% Allgemeine Form:
%    - div(a grad u) + b.grad u + cu = f     in Omega,
%                                  u = dir   auf rand(Omega).
pde = fdpde(grd,opr,fun_a,fun_b,fun_c,fun_f);% Dgl
bcs = fddirichlet(grd,fun_dir);% Randwerte

%% Reduzierte Matrix aufstellen und lösen
uh = fdsolbvp(grd,opr,bcs,pde);

%% Loesung graphisch darstellen
figure(2); clf
   fdplot(grd,uh,gopt);
   title('Numerische Lösung');
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