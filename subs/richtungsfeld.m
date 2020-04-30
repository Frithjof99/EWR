function r = richtungsfeld(bsp)
% CALL:  r = richtungsfeld(bsp);
% EINGABE:
%    bsp .. STRUCT; Beispiel als Struktur
% AUSGABE:
%    r ... STRUCT; Beispiel als Struktur
% BESCHREIBUNG:
% M-file zur Erzeugung eines Richtungsfeldes.
% Fuer die eingegebene Struktur wird das Richtungfeld der Differential-
% gleichung berechnet und zu einem gegebenen Anfangswert die (numerische)
% Loesung eingezeichnet. Ohne Argument wird eine Defaultstruktur
% zurueckgeliefert.
% Die Attribute der Struktur bsp sind:
%    func       Funktionspointer (t,y) \in IRxIR |-> f(t,y) \in IR
%    tmin,tmax  Zeitintervall
%    ymin,ymax  Intervall im Bildraum
%    y0         Startwert zur Zeit tmin
%    N          Aufloesung des Richtungfeldes

% W. Doerfler: Numerische Methoden Dgln, WS 2016/17.
% V 18.10.2016, 26.04.2020.
% Karlsruher Institut fuer Technologie, Fakultaet fuer Mathematik.
% Autor der Ausgangsversion: Christian Kanzow 2006.

%% Analyse data
switch nargin
case 0% Rueckgabe einer Defaultstruktur
   r = struct('func',@(t,y) t.^2-1,'tmin',-2.5,'tmax',2.5, ...
              'ymin',-3,'ymax',3,'y0',-2.2,'N',20,'scale',1);
   return;
case 1% Entpacken
   f = bsp.func;   % Die rechte Seite der ODE
   tmin = bsp.tmin;% Rechteck [tmin,tmax] x [ymin,ymax]
   tmax = bsp.tmax;
   ymin = bsp.ymin;
   ymax = bsp.ymax;
   y0 = bsp.y0;    % Anfangswert
   N  = bsp.N;     % Aufloesung
   r  = bsp;       % Return
   sc = bsp.scale  % Skalierung der Richtungen
otherwise
   error('*** Falsche Anzahl von Parametern ***');
end
fprintf('\nProgramm richtungsfeld\n');tic

%% Sammle Daten
%%% Variante 1: Explizite Schleife
% t = linspace(tmin,tmax,N);
% y = linspace(ymin,ymax,N);
% dy = zeros(N,N);
% dt = zeros(N,N);
% for i=1:length(y)
%    for j=1:length(t)
%       dy(i,j) = f(t(j),y(i));
%       dt(i,j) = 1;
%    end
% end
%%% Variante 2: Implizite Schleife
[t,y] = meshgrid(linspace(tmin,tmax,N),linspace(ymin,ymax,N));
dy = f(t,y);
dt = ones(size(dy));

%% Richtungsbild zeichnen
quiver(t,y,dt,dy,sc);
axis([tmin tmax ymin ymax]);
xlabel('t','FontSize',16); ylabel('y(t)','FontSize',16);
fstr = func2str(f); fstr(1:6)='';
title(['Richtungsfeld fuer die Dgl y'' = ' fstr],'FontSize',14);

%% Loesung berechnen und einzeichnen
hold on
[t,y] = ode23(f,[tmin,tmax],y0);
plot(t,y,'k-','LineWidth',2);
hold off;

%% Ende
fprintf('Fertig\n');toc
