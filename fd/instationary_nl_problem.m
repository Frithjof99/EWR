%% Zusammenfassung: Instationaeres nichtlineares Problem
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
%    d/dt u - div(a grad u) + b.grad u + cu = nlf(u)   in [0,T] x Omega,
%                                         u = dir      auf [0,T] x rand(Omega).
pde = fdpde(grd,opr,fun_a,fun_b,fun_c,fun_f);
bcs = fddirichlet(grd,fun_u0);

%% Reduzierte Matrix aufstellen
L = fdassempde(grd,opr,pde);
[E,ur0] = fdassembcs(grd,opr,bcs);
A = E'*L*E;

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
errLi = zeros(MM,1);% Speichern L^{oo}-Fehler
globt = zeros(MM,1);% Speichern globale Zeit
stp = 4;% Ausgabe mod(.,stp)
ct = 1;% Zaehler
if ~isempty(fun_uext)% Fehlerberechnung
   globt(1) = tbeg;
   errLi(1) = fdnodalerror(grd,E*y0+ur0,@(xx)fun_uext(xx,tbeg));
end
nopt = newton(); nopt.Display = 0; nopt.StepCtrl = 0;

%% Timeloop
yh = E*y0+ur0;% Anfangswert
figure(2); clf;
   fdplot(grd,yh,gopt);
   title(['Loesung zur Zeit t = ' sprintf('%7.4f',tbeg)]);
drawnow
for ii=1:MM
   t = tbeg+ii*dt;
   ur1 = [fun_dir(grd.x(1),t);zeros(grd.ni,1);fun_dir(grd.x(end),t)];
   if isempty(fun_nlf)
      if theta<1.e-8
         y1 = (y0-dt.*(A*y0+E'*L*ur0))-E'*(ur1-ur0);
      else
         y1 = (Id+theta*dt*A)\( y0-E'*(ur1-ur0) ...
                 -dt*( (1-theta)*(A*y0+E'*L*ur0)+theta*E'*L*ur1 ) );
      end
   else
      if theta<1.e-8
         y1 = (y0-dt.*(A*y0+E'*L*ur0-fun_nlf(y0)))-E'*(ur1-ur0);
      else
         NlFunc  = @(vv) vv-(y0-E'*(ur1-ur0) ...
                        -dt*((1-theta)*(A*y0+E'*L*ur0)+theta*(A*vv+E'*L*ur1) ...
                        -fun_nlf((1-theta)*(y0+E'*ur0)+theta*(vv+E'*ur1))));
         NlFuncp = @(vv) Id + dt*theta* ...
            ( A - spdiags(fun_nlfp((1-theta)*(y0+E'*ur0)+theta*(vv+E'*ur1)),0:0,grd.ni,grd.ni) );
         % Newton
         [y1,fu] = newton(NlFunc,NlFuncp,y0,nopt);
      end
   end
   if mod(ii,stp)==0% Ausgabe alle stp Schritte
      yh = E*y1+ur1;% Globale Loesung
      fdplot(grd,yh,gopt);
      if ~isempty(fun_uext)% Fehlerberechnung, Ausgabe
         ct = ct+1;
         globt(ct) = t;% Globale Zeit speichern
         fun_uex = @(xx) fun_uext(xx,t);
         errLi(ct) = fdnodalerror(grd,yh,fun_uex);% Fehler speichern
         if ~isempty(fun_uext)
            hold on; fdplot(grd,fdinterp(grd,opr,fun_uex),gopt2); hold off;
         end
         title(['Loesung zur Zeit t = ' sprintf('%7.4f',t)]);
         drawnow
      end
   end
   y0 = y1;
   ur0 = ur1;
end
yh = E*y1+ur1;% Globale Loesung
fdplot(grd,abs(yh),gopt);
if ~isempty(fun_uext)
   fun_uex = @(xx) fun_uext(xx,tend);
   hold on; fdplot(grd,fdinterp(grd,opr,fun_uex),gopt2); hold off;
end
title(['Loesung zur Zeit t = ' sprintf('%7.4f',ii*dt)]);
drawnow

%% Ausgabe
% Fehler ueber die Zeit
if ~isempty(fun_uext)
   errLi = errLi(1:ct);
   globt = globt(1:ct);
   figure(3); clf
      plot(globt,errLi,'k-','LineWidth',2);
      title('Fehler ueber die Zeit t','FontSize',12);
      xlabel('t','FontSize',12); ylabel('err(t)','FontSize',12);
   drawnow
   errli = max(errLi);
end

%% ENDE
