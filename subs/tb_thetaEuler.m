function [t,y,dtn] = tb_thetaEuler(f,tspan,y0,opts)
%% CALL: [t,y,dtn] = tb_thetaEuler(f,tspan,y0,opts);
% EINGABE:
%    f .. FUNC; f in y'=f(t,y), f:IRxIR^m->IR^m
%    tspan ... DOUBLE(1,2); [tmin,tmax]
%    y0 ... DOUBLE(*); Startwert
%    opts ... STRUCT; Diverse Optionen (s.u.)
% AUSGABE:
%    [t,y] ... DOUBLE(*),DOUBLE(*,m); Vektoren [t_i], [y^i] mit 
%                 y^i\approx y(t_i)
%    dtn ... DOUBLE; Ggf neu berechnetes dt
% BESCHREIBUNG:
% Loesung einer Differentialgleichung mit dem theta-Euler-Verfahren
%   y^{i+1} = y^i
%      + \tau ( (1-\theta) f(t_i,y^i) + \theta f(t_{i+1},y^{i+1}) ).
% Die auftretende nichtlineare Gleichung wird im skalaren Fall mit fzero, im
% Falle eines Systems mit fminsearch oder fsolve geloest. Ohne Parameterliste
% erfolgt die Rueckgabe einer Defaultstruktur fuer opts.

% W. Doerfler: Einfuehrung in das Wissenschaftliche Rechnen, SS 2017. 
% V 25.04.2020.
% Karlsruher Institut fuer Technologie, Fakultaet fuer Mathematik

%% Analysiere Daten
switch nargin
case 0
   optsi = optimset('TolX',1.0e-8,'TolFun',1.0e-10,'MaxFunEvals',2000, ...
                    'MaxIter',2000);
   t = struct('dt',0.1,'theta',1,'ssolver','fsolve','optsi',optsi);
   y = NaN; dtn = NaN;
   return
case 3
   dt = (tspan(2)-tspan(1))/20;
   theta = 1;
   ssolver = 'fsolve';% Systemsolver: fminsearch, fsolve (schneller)
   optsi = optimset('TolX',1.0e-8,'TolFun',1.0e-10,'MaxFunEvals',2000, ...
                    'MaxIter',2000); 
case 4
   dt = opts.dt;
   theta = opts.theta;
   ssolver = opts.ssolver;
   optsi = opts.optsi;
otherwise
end

%% Initialisieren
NN = ceil((tspan(2)-tspan(1))/dt);% Zahl der uniformen Zeitschritte
dt = (tspan(2)-tspan(1))/NN;      % Nachjustieren dt
m = length(y0);                   % Vektordimension
t = zeros(1,NN);
y = zeros(m,NN);
dtn = dt;% Ausgabe
% Festeinstellung
epo = 0;% Parameter (in [0,1]) fuer extrapolierte Startwerte des Nl Loesers

%% Schleife
t(1) = tspan(1);
y(:,1) = y0;
if theta<10^(-10)% Expliziter Fall
   for i=1:NN
      t(i+1) = t(1)+i*dt;
      y(:,i+1) = y(:,i)+dt*f(t(i),y(:,i));
   end
else        % Impliziter Fall
   if m==1  % Skalares y
      for i=1:NN
         t(i+1) = t(1)+i*dt;
         implf = @(z) z-theta*dt*f(t(i+1),z)-y(i) ...
                       -(1-theta)*dt*f(t(i),y(i));
         if epo>10^(-10)
            yextra = y(i)+epo*dt*f(t(i),y(i)); % Extrapolierter Startwert
            %[y(i+1),fval,exfl,out] ...
            y(i+1) = fzero(implf,yextra,optsi);% MATLAB fzero
         else
            y(i+1) = fzero(implf,y(i),optsi);  % MATLAB fzero
         end
      end
   else     % Vektorielles y
      % Optionen sortieren
      switch ssolver
      case 'fsolve'% Feldnamen sind versionsabhaengig!
         optsfz = optimoptions('fsolve');
         optsfz.StepTolerance = optsi.TolX;
         optsfz.FunctionTolerance = optsi.TolFun;
         optsfz.MaxFunctionEvaluations = optsi.MaxFunEvals;
         optsfz.MaxIterations = optsi.MaxIter;
         optsfz.Display = 'none';
         systemsolver = @(fun,y0) fsolve(fun,y0,optsfz);
      case 'fminsearch'
         optsfm = optimset('fminsearch');
         optsfm.TolX = optsi.TolX; optsfm.TolFun = optsi.TolFun;
         optsfm.MaxFunEvals = optsi.MaxFunEvals; optsfm.MaxIter = optsi.MaxIter;
         optsfm.Display = 'none';
         systemsolver = @(fun,y0) fminsearch(@(z)norm(fun(z),2)^2,y0,optsfm);
      end
      for i=1:NN
         t(i+1) = t(1)+i*dt;
         implf = @(z) z-theta*dt*f(t(i+1),z)-y(:,i) ...
                       -(1-theta)*dt*f(t(i),y(:,i));
         if epo>10^(-10)
            yextra = y(:,i)+epo*dt*f(t(i),y(:,i));% Extrapolierter Startwert
            y(:,i+1) = systemsolver(implf,yextra);% Siehe Definition oben
         else
            y(:,i+1) = systemsolver(implf,y(:,i));% Siehe Definition oben
         end
       end
   end
end
%% ENDE