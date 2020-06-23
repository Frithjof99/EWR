%% Zusammenfassung: Eigenwertproblem
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
%    - div(a grad u) + b.grad u + cu = \lambda m u   in Omega,
%                                  u = 0             auf rand(Omega).
pde = fdpde(grd,opr,fun_a,fun_b,fun_c);
bcs = fddirichlet(grd,fun_dir);
pdm = fdpde(grd,opr,[],[],fun_m);

%% Reduzierte Matrizen aufstellen
[A,b,E] = fdassembvp(grd,opr,bcs,pde);
I = fdassembvp(grd,opr,bcs,pdm);

%% Eigenwerte berechnen
% Initialisierung
mm = 1;% mm=0: eig, mm=1: eigs
if size(A,1)<100, mm = 0; end
opts.disp = 0;% Ausgabe?
sigma = 'sm';% Kleinste Eigenwerte
% Eigenwerproblem loesen
if mm==1
   [V,D] = eigs(A,I,nb,sigma,opts);
else
   [V,D] = eig(full(A),full(I));
end
% Eigenwerte fallend sortieren
d = diag(D);
[d,ix] = sort(d);
V = V(:,ix);
% Ueberfluessige Werte abschneiden
if mm==0
   d = d(1:nb);
   V = V(:,1:nb);
end
% Essentielle Randbedingung einarbeiten
U = E*V;

%% Loesung graphisch darstellen
figure(2); clf; hold on% Eigenwerte 1...nb anzeigen
   plot(1:nb,d/pi^2,'ko-','LineWidth',2);
   plot(1,d(1)/pi^2,'ro','LineWidth',2);
   title('Eigenwerte');
   xlabel('j'); ylabel('\lambda_j/\pi^2');
drawnow
figure(3); clf% Eigenfunktionen nacheinander zeigen
   fdplot(grd,real(U(:,1)));
   title('Eigenfunktionen');
   zlabel('u(x,y)');
drawnow
for j=2:nb
   pause(0.5);
   figure(2);% Gezeigte Eigenfunktion markieren
      plot(j-1,d(j-1)/pi^2,'ko','LineWidth',2);
      plot(j,d(j)/pi^2,'ro','LineWidth',2);
      title('Eigenwerte');
      xlabel('j'); ylabel('\lambda_j/\pi^2');
   drawnow
   figure(3);
      fdplot(grd,real(U(:,j)));
      title(['Eigenfunktion ' num2str(j)]);
      zlabel('u(x,y)');
   drawnow
end
figure(2);% Eigenfunktion ent-markieren
   pause(2);
   plot(nb,d(nb)/pi^2,'ko','LineWidth',2);
   hold off
drawnow
figure(3);

%% Ggf. Fehler berechnen
if ~isempty(lam_ex)
   ml = min(length(d),length(lam_ex));
   errlam = max(abs(lam_ex(1:ml)-d(1:ml)))/max(abs(lam_ex(1:ml)));
   fprintf(' N= %d, err1= %7.4e\n',NperDim,errlam);
   errli =  max(abs(lam_ex(1:ml)-d(1:ml)));
end

%% ENDE
