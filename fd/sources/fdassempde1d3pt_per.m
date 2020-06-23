function [L,f] = fdassempde1d3pt_per(grd,pde)
%% CALL: [A,F] = fdassempde1d3pt_per(grd,pde);
% INPUT:
%    grd ... STRUCT; Gitter-Struktur.
%    pde ... STRUCT; Koeffizientenfunktionen [x_1,...,x_d]->IR.
% OUTPUT:
%    L ... DOUBLE**; 3-Punkte-Stern-Matrix.
%    f ... DOUBLE*; Vektor der rechten Seite.
% DESCRIPTION:
% FDASSEMPDE1D3PT_PER Berechnet den 3-Punkte-Stern
%    L = a/h^2 [-1 2 -1] + b/(2h) [-1 0 1] + c [0 1 0]
% und den Lastvektor f auf uniformen 1D Finite-Differenzen-Gittern unter
% Beruecksichtigung periodischer Randbedingungen.
% [A,F] = FDASSEMPDE1D3PT_PER(GRD,OPR,PDE) gibt die Matrix A und den Vektor F der
% rechten Seite des vollen linearen Gleichungssystems zurueck. Das Gitter ist
% durch die Struktur GRD mit den Komponenten G, X, ... definiert (siehe FDGRID),
% die Funktionspointer der Koeffizientenfunktionen und der rechten Seite der
% partiellen Dgl befinden sich in der Struktur PDE.
% FDASSEMPDE1D3PT_PER wird direkt aufgerufen.

% Version 1.0: Willy Doerfler, KIT, 2013.

%% Initialisierungen
G = grd.G;
x = grd.x;
N = grd.nn;% Die Anzahl aller Gitterpunkte
h = x(2)-x(1);% Bestimme die Gitterweite (Annahme: uniformes Gitter)
h2i = 1/h^2;% Diagonalskalierung
%
f = zeros(N,1);% Allokiere Speicher fuer die rechte Seite
% Vektoren, in denen die Zeilen- und Spaltenindizes sowie die Werte der
% Systemmatrix gespeichert werden
k = zeros(1,3*N);
l = zeros(1,3*N);
s = zeros(1,3*N);
ffun = pde.loadfuns;
afun = pde.stifcoefs;
bfun = pde.trnscoefs;
cfun = pde.masscoefs;
% Faelle: func=[] oder func=const zu func(xx)
if isempty(ffun),       ffun = @(xx) 0;
elseif isnumeric(ffun), ffun = @(xx) ffun; end
if isempty(afun),       afun = @(xx) 0;
elseif isnumeric(afun), afun = @(xx) afun; end
if isempty(bfun),       bfun = @(xx) 0;
elseif isnumeric(bfun), bfun = @(xx) bfun; end
if isempty(cfun),       cfun = @(xx) 0;
elseif isnumeric(cfun), cfun = @(xx) cfun; end

%% Wende den 3-pt-Stern an
ic = 0;
for i=1:N
   gi = G(i);
   if gi==0, continue; end% Naechster Index
   if i<N                 % Innerer Punkt
      % Rechte Seite
      f(gi) = ffun(x(i));

      % Koeffizienten
      aval = afun(x(i));
      bval = bfun(x(i));
      cval = cfun(x(i));
      bh = bval*h/2;% Zentriert fuer b

      % Zentrum des 3-Punkte-Sterns
      ic = ic + 1;
      k(ic) = gi;
      l(ic) = gi;
      s(ic) = 2*aval*h2i+cval;

      % Linker Ast des 3-Punkte Sterns 
      ic = ic + 1;
      k(ic) = gi;
      if i==1
         l(ic) = G(N-1);
      else
         l(ic) = G(i-1);
      end
      s(ic) = h2i*(-aval-bh);
            
      % Rechter Ast des 3-Punkte-Sterns
      ic = ic + 1;
      k(ic) = gi;
      if i==N-1
         l(ic) = G(1);
      else
         l(ic) = G(i+1);
      end
      s(ic) = h2i*(-aval+bh);
   else                   % Linker Randpunkt
      ic = ic + 1;
      k(ic) = gi;
      l(ic) = gi;
      s(ic) = 1;
   end
end

%% Erzeuge die Systemmatrix
L = sparse(k(1:ic),l(1:ic),s(1:ic),N,N);
% Strip
L(:,end) = [];
L(end,:) = [];
f(end) = [];

return
