function opt = fdplot(grd,vh,opt)
%% CALL: opt = fdplot(grd,vh,opt);
% INPUT:
%    grd ... STRUCT; Gitter-Struktur (fdgrid).
%    vh ... DOUBLE(*,d); Gitterfunktion [DEF 0], d=2,3: Vektorfeld.
%    opt ... STRUCT; Options-Struktur [DEF by return].
% OUTPUT:
%    opt ... STRUCT; Options-Struktur mit fdplot-defaults fuer nargin==0.
% DESCRIPTION:
% FDPLOT Graphische Ausgabe des Rechengebietes oder einer Gitterfunktion.
% FDPLOT() gibt eine Struktur mit den Defaultoptionen zurueck. Mit
%    diesen lassen sich Label, Linientyp, Ausgabeart oder eine feste Achse
%    waehlen.
% FDPLOT(GRD) veranlasst die graphische Ausgabe des Rechengitters:
%    Innere Punkte sind blau, Randpunkte rot dargestellt.
% FDPLOT(GRD,VH) veranlasst die graphische Ausgabe der Gitterfunktion VH mit
%    Defaultoptionen.
% FDPLOT(GRD,VH,OPT) verwendet die uebergebenen Optionen.
% Struktur OPT:
%    XLBL, YLBL, ZLBL ... CHAR; Benennung Label [DEF 'x' 'y' 'z']
%    LTYP ... CHAR; Linientyp [DEF 'k-']
%    LWID ... INT; Linienbreite [DEF 2]
%    FIXAXIS ... Feste Achsen [DEF []]
%    METH2D ... Grafikmethode in 2D [DEF surf] {surfc, contour},
%           ... vh(*,d), d=2,3: Quiver.
%    SFUN ... FUNC; Zeichne SFUN(VH) [DEF @(z) z]
%%

% Version 1.0: Willy Doerfler, KIT, 2014/15.
% Version 2.0: Willy Doerfler, KIT, 2020.

%% Unterscheide Eingabedaten
if nargin==0
   opt = defaultplotoption();
   return
elseif nargin==1
   opt = defaultplotoption();
   if grd.dim==1
      x = grd.x;
      LL = x(end)-x(1);
      hold on;
      plot(x,zeros(size(x)),'b.');
      plot([x(1),x(end)],[0 0],'r.');
      axis equal;
      xlabel(opt.xlbl,'FontSize',12);
      ylabel('');
      axis([x(1)-0.05*LL x(end)+0.05*LL -0.1*LL 0.1*LL]);
      set(gca,'yticklabel',{''});
      hold off;
   elseif grd.dim==2
      [i,j] = find(grd.G);% Alle Knoten des Gebietes
      plot(grd.x(i),grd.y(j),'marker','*','linestyle','none','color','b');
      i = grd.bverts(1,:);% Randknoten
      j = grd.bverts(2,:);
      hold on;
      plot(grd.x(i),grd.y(j),'marker','*','linestyle','none','color','r');
      xlabel(opt.xlbl,'FontSize',12);
      ylabel(opt.ylbl,'FontSize',12);
      hold off;
   elseif grd.dim==3
      fdplot3(grd);
   end
elseif nargin>=2
   if nargin==2, opt = defaultplotoption(); end
   if nargin==3 && isempty(opt), opt = defaultplotoption(); end
   if grd.dim==1
      plot(grd.x,opt.sfun(vh(grd.G)),opt.ltyp,'LineWidth',opt.lwid);
      if ~isempty(opt.fixaxis), axis(opt.fixaxis(1:4)); end
      xlabel(opt.xlbl,'FontSize',12);
      ylabel(opt.ylbl,'FontSize',12);
   elseif grd.dim==2
      if size(vh,2)==1% Skalar
         [X,Y,Z] = fdmeshgrid(grd,opt.sfun(vh));
         if strcmp(opt.meth2d,'surfc')
            surfc(X,Y,Z);
         elseif strcmp(opt.meth2d,'contour')
            contour(X,Y,Z);
         else
            surf(X,Y,Z);
         end
      elseif size(vh,2)==2% Vektorfeld IR^2
         [X,Y,Zall] = fdmeshgrid(grd,opt.sfun(vh(:,1:2)));
         quiver(X,Y,Zall(:,:,1),Zall(:,:,2));
      elseif size(vh,2)==3% Vektorfeld IR^3
         [X,Y,Zall] = fdmeshgrid(grd,opt.sfun(vh(:,1:3)));
         Z = sparse(size(X,1),size(X,2));
         quiver3(X,Y,Z,Zall(:,:,1),Zall(:,:,2),Zall(:,:,3), ...
                 opt.qscal,'LineWidth',opt.qlw,'MarkerSize',opt.qms);
      end
      if ~isempty(opt.fixaxis), axis(opt.fixaxis); end
      if ~isempty(opt.view), view(opt.view); end
      xlabel(opt.xlbl,'FontSize',12);
      ylabel(opt.ylbl,'FontSize',12);
      zlabel(opt.zlbl,'FontSize',12);
   elseif grd.dim==3
      fdplot3(grd,vh,opt);
   end
end
if ~isempty(opt.titlestr), title('opt.titlestr'); end
opt = [];% Return []

return

%% ------------------------------------------------------------------------

function opt = defaultplotoption()
   opt = struct('xlbl','x','ylbl','y','zlbl','z','ltyp','k-','lwid',2, ...
            'titlestr',[],'fixaxis',[],'view',[],'meth2d','surf', ...
            'sfun',@(z) z,'qscal',2,'qlw',1.2,'qms',0.1, ...
            'xslice',0,'yslice',0,'zslice',0);
return
