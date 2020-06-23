function [L,f] = fdassempde2d9pt(grd,pde)
%% CALL: [L,f] = fdassempde2d9pt(grd,pde);
% INPUT:
%    grd ... STRUCT; Gitter-Struktur.
%    pde ... STRUCT; Koeffizientenfunktionen.
% OUTPUT:
%    L ... DOUBLE**; 9-Punkte-Stern-Matrix.
%    f ... DOUBLE*; Vektor der rechten Seite.
% DESCRIPTION:
% FDASSEMPDE2D9PT berechnet den 9-Punkte-Stern
%        [-1/6  -4/6    -1/6 ]
%    L = [-4/6  (20/6)  -4/6 ] *a/h^2, a=konstant
%        [-1/6  -4/6    -1/6 ]
% und den Lastvektor f auf uniformen 2D Finite-Differenzen-Gittern.
% [A,F] = FDASSEMPDE2D9PT(GRD,OPR,PDE) gibt die Matrix A und den Vektor F der
% rechten Seite des vollen linearen Gleichungssystems zurueck. Das Gitter ist
% durch die Struktur GRD mit den Komponenten G, X, ... definiert (siehe FDGRID),
% die Funktionspointer der Koeffizientenfunktionen und der rechten Seite der
% partiellen Dgl befinden sich in der Struktur PDE.
% FDASSEMPDE2D9PT wird nicht direkt, sondern nur ueber FDASSEMPDE aufgerufen.

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
afun = pde.stifcoefs;
if isnumeric(afun)
   a = afun;
else
   error(' *** a sollte konstant sein ***');
end
Lii =  20*a/6/h^2;% Diagonalelement
Lij = -4*a/6/h^2;% Nebendiagonalelement
Liijj = -a/6/h^2;% Eckenelement
%
f = zeros(N,1);% Allokiere Speicher fuer die rechte Seite
% Vektoren, in denen die Zeilen- und Spaltenindizes sowie die Werte der
% Systemmatrix gespeichert werden
k = zeros(1,9*N);
l = zeros(1,9*N);
s = zeros(1,9*N);
rhsfun = pde.loadfuns;

%% Markiere Randknoten
for ir=1:NB
   i = grd.bverts(1,ir);
   j = grd.bverts(2,ir);
   G(i,j) = -G(i,j);
end

%% Wende den 9-pt-Stern an
ic = 0;
for i=1:m
   for j=1:n
      if G(i,j)==0, continue; end% Naechster Index
      if G(i,j)>0                % Innerer Punkt
         % Rechte Seite
         f(G(i,j)) = ( 8*rhsfun(x(i),y(j)) ...
            +rhsfun(x(i-1),y(j))+rhsfun(x(i+1),y(j)) ...
            +rhsfun(x(i),y(j-1))+rhsfun(x(i),y(j+1)) )/12;

         % Zentrum des 9-Punkte-Sterns
         ic = ic + 1;
         k(ic) = G(i,j);
         l(ic) = G(i,j);
         s(ic) = Lii;

         % Linker Ast des 9-Punkte Sterns 
         ic = ic + 1;
         k(ic) = G(i,j); 
         l(ic) = G(i-1,j);
         s(ic) = Lij;
            
         % Rechter Ast des 9-Punkte-Sterns
         ic = ic + 1;
         k(ic) = G(i,j);
         l(ic) = G(i+1,j);
         s(ic) = Lij;

         % Unterer Ast des 9-Punkte-Sterns
         ic = ic + 1;
         k(ic) = G(i,j);
         l(ic) = G(i,j-1);
         s(ic) = Lij;

         % Oberer Ast des 9-Punkte-Sterns
         ic = ic + 1;
         k(ic) = G(i,j);
         l(ic) = G(i,j+1);
         s(ic) = Lij;

         % Linke obere Ecke des 9-Punkte-Sterns
         ic = ic + 1;
         k(ic) = G(i,j);
         l(ic) = G(i-1,j+1);
         s(ic) = Liijj;

         % Rechte obere Ecke des 9-Punkte-Sterns
         ic = ic + 1;
         k(ic) = G(i,j);
         l(ic) = G(i+1,j+1);
         s(ic) = Liijj;

         % Linke untere Ecke des 9-Punkte-Sterns
         ic = ic + 1;
         k(ic) = G(i,j);
         l(ic) = G(i-1,j-1);
         s(ic) = Liijj;

         % Rechte untere Ecke des 9-Punkte-Sterns
         ic = ic + 1;
         k(ic) = G(i,j);
         l(ic) = G(i+1,j-1);
         s(ic) = Liijj;
      else              % Randpunkt
         ic = ic + 1;
         k(ic) = -G(i,j);
         l(ic) = -G(i,j);
         s(ic) = 1;
      end
   end
end

%% Erzeuge die Systemmatrix
L = sparse(abs(k(1:ic)),abs(l(1:ic)),s(1:ic),N,N);

return
