function grd = fdgrid(domain,m,geom)
%% CALL: grd = fdgrid(domain,m,geom);
% INPUT:
%    domain ... STRING; Symbol des vorimplementierten Gebietes (s.u.)
%               [DEF 'S'].
%    m ... INT; Anzahl der Gitterpunkte je Dimension [DEF 20].
%    geom ... DOUBLE*; Geometrie des Rechengebietes (s.u.)
%             [DEF [zeros(dim) 2*ones(dim)]].
% OUTPUT:
%    grd ... STRUCT; Gitter-Struktur.
% DESCRIPTION:
% FDGRID Finite-Differenzen-Gitter auf Standardrechengebieten.
% GRD = FDGRID(REGION,M,A) gibt eine Struktur mit den Attributen
%    G, X, Y, BVERTS, NN, NI, DIM, ORDERING
% zurück. G>0 gibt das innere Gebiet an und BVERTS listet die Randknoten.
%    BVERTS(IR) = [I,J,BTYP],
% heisst: Randknoten IR hat die Nummer G(I,J) und den Typ BTYP (BTYP=0 etwa
% bedeutet Dirichletknoten, BTYP=1 Neumannknoten). Im Fall des Torus gibt BTYP
% die Randknotennummer des identifizierten Punktes.
% X, Y (1d: X) sind die Koordinatenvektoren des Finite-Differenzen-
% Gitters.
% Das Rechengebiet wird durch den Eingabeparameter DOMAIN bestimmt. 
% Folgende Rechengebiete sind z.B. verfügbar:
%     'I' - Intervall (DIM=1)
%     'S' - Ein Quadrat mit Schwerpunkt im Koordinatenursprung (DIM=2),
%     'L' - Das L-Gebiet mit einspringender Ecke in [0,0] (DIM=2),
%     'D' - Eine Kreisscheibe mit Mittelpunkt in [0,0] (DIM=2),
%     'T' - Torus, 'S' mit Randidentifikation (DIM=2).
%     'W' - Wuerfel mit Schwerpunkt im Koordinatenursprung (DIM=3).
% Weitere Gebiete fuer DIM=2: Siehe NUMGRID (A,B,C,H,N).
% Der Eingabeparameter M legt die Anzahl der Gitterintervalle je Dimension
% fest. Der Eingabeparameter GEOM bestimmt die Lage des Rechengebietes durch
% Angabe von Zentrum und Ausdehnung:
%     [C SIZE] - C=center in IR^d, SIZE=[xdiam,...] in IR^d.

% Version 0.0: Markus Richter, KIT, 2008.
% Version 1.0: Willy Doerfler, KIT, May 2020.

%% Analysiere Eingabedaten
if nargin<3
   geom = [];
end
if nargin<2
   m = 20;
end
if nargin<1
   domain = 'S';
end
if (strcmp(domain,'L') && mod(m,2)==0)% Gerades m fuer das Gebiet 'L'
   m = m+1;
   disp([mfilename ' Increased m by one']);
end
if strcmp(domain,'I'), dim = 1; end
if strcmp(domain,'S') || strcmp(domain,'L') || strcmp(domain,'D') ...
                      || strcmp(domain,'T')
   dim = 2;
end
if strcmp(domain,'W'), dim = 3; end

%% DIM=1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dim==1
   if isempty(geom), geom = [0 2]; end
   G = 1:m;
   bverts = [1 m; 0 0];
   x = linspace(geom(1)-geom(2)/2,geom(1)+geom(2)/2,m);
   grd = struct('G',G,'x',x,'bverts',bverts,'nn',m,'ni',m-2, ...
                'dim',1,'ordering','natural');
   return
end

%% DIM=2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dim==2
%% Grid aus numgrid (innere Punkte)
if strcmp(domain,'T'), domain = 'S'; per = 1; else per = 0; end% Case Torus
if isempty(geom), geom = [0 0 2 2]; end
G = numgrid(domain,m+2);% Etwas groesser waehlen!
if (strcmp(domain,'L'))% L-Gebiet nachbessern (innerer Rand auf y=0)
   G((m+3)/2,2:end-1) = 1; G(2:end-1,(m+3)/2) = 1;
end
G = gridrenumber(G);% Numerierung ist unten-links nach oben rechts
% G_{ij} soll nun transformiert werden, so dass G_{ij} mit [x(i),y(j)]
% korrespondiert!
G = flipud(G)';

%% Randpunkte setzen
nn = max(G(:));% Knotenanzahl
nb = 0;% Init Anzahl Randpunkte
[mg,ng] = size(G);
bverts = zeros(3,nn);
for j=2:ng-1% Schleife ueber alle Punkte (2:m+1)x(2:m+1)
   jma = j-1;%
   jmi = j+1;
   for i=2:mg-1
      ima = i-1;
      imi = i+1;
      if G(i,j)>0% Innerer Punkt ...
         if G(ima,j)==0 || G(imi,j)==0 || ... % ... mit Nachbarn aussen
               G(i,jma)==0 || G(i,jmi)==0 || ...
               G(ima,jma)==0 || G(imi,jma)==0 || ...
               G(ima,jmi)==0 || G(imi,jmi)==0
            nb = nb+1;
            bverts(1,nb) = i;
            bverts(2,nb) = j;
            bverts(3,nb) = 0;
         end
      end
   end
end
ni = nn-nb;% Anzahl innerer Punkte
bverts = bverts(:,1:nb);% Ueberfluessigen Speicher freigeben
G = G(2:mg-1,2:ng-1);% Ueberfluessigen Rand abschneiden
bverts(1:2,:) = bverts(1:2,:)-1;% Entspr. Indexkorrektur in bverts
x = linspace(geom(1)-geom(3)/2,geom(1)+geom(3)/2,m);% Koordinaten
y = linspace(geom(2)-geom(4)/2,geom(2)+geom(4)/2,m);
% Case Torus: Setze Identifikationspunkte: bvert(3,:) ist Randknotennummer
% des identifizierten Punktes.
if per
   N = size(G,1);%
   r0 = 3*N-4;
   for jj=1:N% Id letzte Spalte
      rnr = r0+jj;
      bverts(3,rnr) = jj;
   end
   jj = 1;
   rnr = N+2*(jj-1);
   bverts(3,rnr) = rnr-N+1;
   for jj=2:N-1% Id letzte Zeile
      rnr = N+2*(jj-1);
      bverts(3,rnr) = rnr-1;
   end
   rnr = rnr+N;
   bverts(3,rnr) = 1;
end

%% Struktur belegen
grd = struct('G',G,'x',x,'y',y,'bverts',bverts,'nn',nn,'ni',ni, ...
             'dim',2,'ordering','lexicographic','id','fd');
end

%% DIM=3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dim==3
   %grd = fdgrid3(domain,m,geom,implicit_grid);% ToDo
   grd = fdgrid3(domain,m,geom);
end

return

%% ------------------------------------------------------------------------

function G = gridrenumber(G)
% G von 'numgrid' ist spaltenweise von oben nach unten numeriert.
% Damit es besser zu den Koordinatenfeldern x,y passt, wird G
% umnummeriert: Von links unten nach rechts oben.
[mg,ng] = size(G);
ni = 0;
for i=mg:-1:1
   for j=1:ng
      if G(i,j)>0
         ni = ni+1;
         G(i,j) = ni;
      end
   end
end

return
