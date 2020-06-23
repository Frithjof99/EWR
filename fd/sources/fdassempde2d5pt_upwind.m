function [L,f] = fdassempde2d5pt_upwind(grd,pde)
%% CALL: [L,f] = fdassempde2d5pt_upwind(grd,pde);
% INPUT:
%    grd ... STRUCT; Gitter-Struktur.
%    pde ... STRUCT; Koeffizientenfunktionen [x_1,...,x_d]->IR.
% OUTPUT:
%    L ... DOUBLE**; 5-Punkte-Stern-Matrix.
%    f ... DOUBLE*; Vektor der rechten Seite.
% DESCRIPTION:
% FDASSEMPDE2D5PT_UPWIND berechnet den 5-Punkte-Stern
%        [                      -a+b2m*h                     ]
%    L = [ -a-b1p*h  4*a+(b1p-b1m+b2p-b2m)*h+c*h^2  -a+b1m*h ] /h^2
%        [                      -a-b2p*h                     ]
% (mit b1p=max(b1,0), b1m=min(b1,0), b2p=max(b2,0), b2m=min(b2,0)) und den
% Lastvektor f auf uniformen 2D Finite-Differenzen-Gittern.
% [A,F] = FDASSEMPDE2D5PT_UPWIND(GRD,OPR,PDE) gibt die Matrix A und den Vektor F
% der rechten Seite des vollen linearen Gleichungssystems zurueck. Das Gitter
% ist durch die Struktur GRD mit den Komponenten G, X, ... definiert (siehe
% FDGRID), die Funktionspointer der Koeffizientenfunktionen und der rechten
% Seite der partiellen Dgl befinden sich in der Struktur PDE.
% FDASSEMPDE2D5PT_UPWIND wird nicht direkt, sondern nur ueber FDASSEMPDE
% aufgerufen.

% Version 0.0: Markus Richter, KIT, 2008.
% Version 1.0: Willy Doerfler, KIT, 2013.

%% Initialisierungen
G = grd.G;
x = grd.x;
y = grd.y;
[m,n] = size(G);% Bestimme die Dimensionen des FD-Gitters
N  = grd.nn;% Die Anzahl aller Gitterpunkte
NI = grd.ni;% Die Anzahl der inneren Gitterpunkte
NB = N-NI;% Die Anzahl der Gitterpunkte am Rand
h  = x(2)-x(1);% Bestimme die Gitterweite (Annahme: uniformes Gitter)
h2i = 1/h^2;% Diagonalskalierung
%
f = zeros(N,1);% Allokiere Speicher fuer die rechte Seite
% Vektoren, in denen die Zeilen- und Spaltenindizes sowie die Werte der
% Systemmatrix gespeichert werden
k = zeros(1,5*N);
l = zeros(1,5*N);
s = zeros(1,5*N);
ffun = pde.loadfuns;
afun = pde.stifcoefs;
bfun = pde.trnscoefs;
cfun = pde.masscoefs;
% Faelle: func=[] oder func=const zu func(xx,yy)
if isempty(ffun),       ffun = @(xx,yy) 0;
elseif isnumeric(ffun), ffun = @(xx,yy) ffun; end
if isempty(afun),       afun = @(xx,yy) 0;
elseif isnumeric(afun), afun = @(xx,yy) afun; end
if isempty(bfun),       bfun = @(xx,yy) [0,0];
elseif isnumeric(bfun), bfun = @(xx,yy) bfun; end
if isempty(cfun),       cfun = @(xx,yy) 0;
elseif isnumeric(cfun), cfun = @(xx,yy) cfun; end

%% Markiere Randknoten
for ir=1:NB
   i = grd.bverts(1,ir);
   j = grd.bverts(2,ir);
   G(i,j) = -G(i,j);
end

%% Wende den 5-pt-Stern an
ic = 0;
for i=1:m
   for j=1:n
      if G(i,j)==0, continue; end% Naechster Index
      if G(i,j)>0                % Innerer Punkt
         % Rechte Seite
         f(G(i,j)) = ffun(x(i),y(j));

         % Koeffizienten
         aval = afun(x(i),y(j));
         bval = bfun(x(i),y(j));
         cval = cfun(x(i),y(j));
         bh1p = max(bval(1),0)*h;% Upwind fuer b1+
         bh1m = min(bval(1),0)*h;% Upwind fuer b1-
         bh2p = max(bval(2),0)*h;% Upwind fuer b2+
         bh2m = min(bval(2),0)*h;% Upwind fuer b2-

         % Zentrum des 5-Punkte-Sterns
         ic = ic + 1;
         k(ic) = G(i,j);
         l(ic) = G(i,j);
         s(ic) = h2i*(4*aval+bh1p-bh1m+bh2p-bh2m)+cval;

         % Linker Ast des 5-Punkte Sterns 
         ic = ic + 1;
         k(ic) = G(i,j); 
         l(ic) = G(i-1,j);
         s(ic) = h2i*(-aval-bh1p);
            
         % Rechter Ast des 5-Punkte-Sterns
         ic = ic + 1;
         k(ic) = G(i,j);
         l(ic) = G(i+1,j);
         s(ic) = h2i*(-aval+bh1m);

         % Unterer Ast des 5-Punkte-Sterns
         ic = ic + 1;
         k(ic) = G(i,j);
         l(ic) = G(i,j-1);
         s(ic) = h2i*(-aval-bh2p);

         % Oberer Ast des 5-Punkte-Sterns
         ic = ic + 1;
         k(ic) = G(i,j);
         l(ic) = G(i,j+1);
         s(ic) = h2i*(-aval+bh2m);
      else              % Randpunkt
         ic = ic + 1;
         k(ic) = G(i,j);
         l(ic) = G(i,j);
         s(ic) = 1;
      end
   end
end

%% Erzeuge die Systemmatrix
L = sparse(abs(k(1:ic)),abs(l(1:ic)),s(1:ic),N,N);

end
