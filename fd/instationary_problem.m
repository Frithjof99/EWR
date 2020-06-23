%% Zusammenfassung: Instationaeres Problem
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
%    d/dt u - div(a grad u) + b.grad u + cu = f     in Omega,
%                                         u = dir   auf rand(Omega).
pde = fdpde(grd,opr,fun_a,fun_b,fun_c,fun_f);
bcs = fddirichlet(grd,fun_dir);

%% Reduzierte Matrix aufstellen
[A,b,E,ur] = fdassembvp(grd,opr,bcs,pde);

%% Initialisierungen
h = grd.x(2)-grd.x(1);
if strcmp(dtscal,'cfl')% Setzung der Zeitschrittweite dt sim h
   dt = dtcon*h;
elseif strcmp(dtscal,'diff')% Setzung der Zeitschrittweite dt sim h*h
   dt = dtcon*h*h;
else
   error('*** Zeitschritt setzen ***');
end
MM = ceil((tend-tbeg)/dt);% Zahl der Zeitschritte
dt = (tend-tbeg)/MM;% Nachjustieren dt
y0 = E'*fdinterp(grd,opr,fun_u0);% Anfangsfunktion interpolieren
y1 = zeros(grd.ni,1);
Id = speye(size(A));
stp = 20;% Ausgabe mod(.,stp)

%% Timeloop
yh = E*y0+ur;% Anfangswert
figure(2); clf;
   fdplot(grd,yh,gopt);
   title(['Loesung zur Zeit t = ' sprintf('%7.4f',tbeg)]);
if ~isempty(fun_uext)% Fehlerberechnung
   if grd.dim==1
      fun_uex = @(xx) fun_uext(xx,tbeg);
      hold on; fdplot(grd,fdinterp(grd,opr,fun_uex),gopt2); hold off;
   elseif grd.dim==2
      fun_uex = @(xx,yy) fun_uext(xx,yy,tbeg);
   end
   globt = zeros(MM+1,1);% Speichern globale Zeit
   errLi = zeros(MM+1,1);% Speichern L^{oo}-Fehler
   ct = 1;% Zaehler
   globt(ct) = tbeg;% Globale Zeit speichern
   errLi(ct) = fdnodalerror(grd,yh,fun_uex);% Fehler speichern
end
drawnow
for ii=1:MM
   t = tbeg+ii*dt;
   if theta<1.e-8% Expliziter Euler
      y1 = y0-dt*(A*y0-b);
   else          % Implizite Euler-Varianten
      y1 = (Id+theta*dt*A)\(y0-dt*((1-theta)*A*y0-b));
   end
   if mod(ii,stp)==0% Ausgabe alle stp Schritte
      yh = E*y1+ur;% Globale Loesung
      fdplot(grd,yh,gopt);
      title(['Loesung zur Zeit t = ' sprintf('%7.4f',t)]);
      if ~isempty(fun_uext)% Fehlerberechnung, Ausgabe
         if grd.dim==1
            fun_uex = @(xx) fun_uext(xx,t);
            hold on; fdplot(grd,fdinterp(grd,opr,fun_uex),gopt2); hold off;
         elseif grd.dim==2
            fun_uex = @(xx,yy) fun_uext(xx,yy,t);
         end
         ct = ct+1;
         globt(ct) = t;% Globale Zeit speichern
         errLi(ct) = fdnodalerror(grd,yh,fun_uex);% Fehler speichern
      end
      drawnow
   end
   y0 = y1;
end
yh = E*y1+ur;% Globale Loesung
fdplot(grd,yh,gopt);
   title(['Loesung zur Zeit t = ' sprintf('%7.4f',ii*dt)]);
drawnow

%% Ausgabe
% Fehler ueber die Zeit
if ~isempty(fun_uext)
   if grd.dim==1
      fun_uex = @(xx) fun_uext(xx,tend);
      hold on; fdplot(grd,fdinterp(grd,opr,fun_uex),gopt2); hold off;
   elseif grd.dim==2
      fun_uex = @(xx,yy) fun_uext(xx,yy,tend);
   end
   if globt(ct)<tend% Letzten Zeitpunkt dazufuegen
      ct = ct+1;
      globt(ct) = tend;% Globale Zeit speichern
      errLi(ct) = fdnodalerror(grd,yh,fun_uex);% Fehler speichern
   end
   errLi = errLi(1:ct);
   globt = globt(1:ct);
   figure(3); clf
      plot(globt,errLi,'k-','LineWidth',2);
      title('Nodaler Fehler ueber die Zeit t','FontSize',12);
      xlabel('t','FontSize',12); ylabel('err(t)','FontSize',12);
   drawnow
   errli = max(errLi);
end

%% ENDE

